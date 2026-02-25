library(Matrix)
library(dplyr)

transcript_meta_analysis <- function(
  study_path_vector,             # character vector of .RData files, each defining `assoc`
  grouping_path_vector,          # character vector of transcript grouping .RData files
  test            = c("Burden", "SKAT", "SKATO", "SMMAT"),
  min_study_cmac  = 1,           # per-study MAC threshold for including a transcript
  min_meta_cmac   = 20,          # min MAC for Cauchy meta p-value combination
  rho             = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1),
  use.anytranscript = TRUE,      # if TRUE, also treat each group_id as its own transcript
  combine.pval      = TRUE,      # if TRUE, call comb_meta_pvlaue at the end
  progress_bar      = TRUE       # show a progress bar over genes
) {

  ## 1. Load per-study association objects (`assoc`) ----------------------
  n_studies  <- length(study_path_vector)
  study_list <- vector("list", n_studies)

  for (i in seq_len(n_studies)) {
    load(study_path_vector[i])          # must create object named `assoc`
    study_list[[i]] <- assoc
  }

  ## 2. Collect all genes (groupings) across studies ----------------------
  groupings <- NULL
  for (i in seq_len(n_studies)) {
    groupings <- unique(c(groupings, study_list[[i]]$results$genename))
  }

  ## 3. Load transcript grouping info for each study ----------------------
  transgroup_list <- vector("list", n_studies)

  for (i in seq_len(n_studies)) {
    if (use.anytranscript) {
      grp_name <- load(grouping_path_vector[[i]])
      grp0     <- get(grp_name)     # original grouping
      grp1     <- grp0              # copy to create "any transcript" rows

      # Normalize the transcript column name to `transcript_id`
      names(grp1)[names(grp1) %in% "TranscriptID"] <-
        names(grp0)[names(grp0) %in% "TranscriptID"] <- "transcript_id"

      # Any-transcript: treat group_id itself as transcript_id
      grp1$transcript_id <- grp1$group_id

      # Stack original and "any transcript" mapping
      transgroup_list[[i]] <- rbind(grp0, grp1)
    } else {
      grp_name <- load(grouping_path_vector[[i]])
      transgroup_list[[i]] <- get(grp_name)
    }
  }

  ## 4. Loop over genes (groupings) --------------------------------------
  res <- NULL
  num <- 1L
  n__groupings <- length(groupings)

  # --- progress bar over genes ---
  pb <- NULL
  if (progress_bar) {
    pb <- utils::txtProgressBar(min = 0, max = n__groupings, style = 3)
  }

  for (group in groupings) {

    # update progress bar
    if (progress_bar && !is.null(pb)) {
      utils::setTxtProgressBar(pb, num)
    }

    # Simple progress messages every 10 genes (and at start/end)
    outnums <- c(1, seq(0, n__groupings, by = 10), n__groupings)
    #cat(num, "\n")
    #if (num %in% outnums) {
    #  cat(num, "out of", n__groupings, "is running\n")
    #}
    num <- num + 1L

    ## 4a. Collect all transcripts for this gene across studies ----------
    transcripts <- NULL
    for (i in seq_len(n_studies)) {
      genegroup  <- subset(transgroup_list[[i]], group_id == group)
      transcripts <- unique(c(transcripts, genegroup$transcript_id))
      transcripts <- transcripts[!is.na(transcripts)]
    }

    if (length(transcripts) == 0L) next

    ## 5. Loop over transcripts for this gene ----------------------------
    for (tt in seq_along(transcripts)) {
      #cat("  transcript index:", tt, "\n")

      transID <- transcripts[tt]

      # Accumulators for this transcript across studies
      variant.id <- NULL     # all variant IDs contributing to this transcript
      n.alts     <- NULL     # per-study total MAC
      V.list     <- vector("list", n_studies)   # per-study covariance
      sv.list    <- vector("list", n_studies)   # per-study score table

      ## 5a. Per-study extraction and SPA-calibrated covariance ----------
      for (i in seq_len(n_studies)) {

        transgroup_data      <- transgroup_list[[i]]
        trans_group_variants <- subset(transgroup_data, transcript_id == transID)

        if (nrow(trans_group_variants) == 0L) next

        # Build variant.id in the same way as in variantInfo
        trans_group_variants$variant.id <- paste(
          trans_group_variants$chr,
          trans_group_variants$pos,
          trans_group_variants$ref,
          trans_group_variants$alt,
          sep = ":"
        )

        total_variant <- study_list[[i]]$variantInfo[[group]]
        total_cov     <- study_list[[i]]$covariance_matrix[[group]]

        # Subset per-study variant info to those in this transcript
        sub_score <- total_variant[total_variant$variant.id %in%
                                     trans_group_variants$variant.id, ]
        sub_score <- unique(sub_score)

        # Skip if no variants or MAC below threshold
        if (is.null(sub_score)) next 
        if (nrow(sub_score) == 0L) next
        if (sum(sub_score$MAC) < min_study_cmac) next

        # Track variant IDs and MAC across studies
        variant.id <- c(variant.id, sub_score$variant.id)
        n.alts     <- c(n.alts, sum(sub_score$MAC))

        # Order by variant.id
        sub_score <- sub_score[order(sub_score$variant.id), ]
        vindex    <- sub_score$variant.id

        # Store per-study score table
        sv.list[[i]] <- sub_score

        # Extract study covariance for these variants
        V.list[[i]] <- total_cov[vindex, vindex, drop = FALSE]

        # Replace diagonal with SPA-calibrated variances (keep off-diagonals)
        if ("SPA.Score.Variance" %in% names(sub_score)) {

        if (nrow(sub_score) > 1L) {
          diag(V.list[[i]]) <- sub_score$SPA.Score.Variance
        } else {
          V.list[[i]][1, 1] <- sub_score$SPA.Score.Variance
        }
        }
      } # end per-study loop

      # Build a summary row for this transcript (n.studies.contributing
      # will be corrected after identifying effective_studies)
      out <- data.frame(
        Group                  = group,
        Transcript             = transID,
        n.studies.contributing = length(sv.list),         # placeholder
        n.site                 = length(unique(variant.id)),
        n.alt                  = sum(n.alts),
        stringsAsFactors       = FALSE
      )
      rownames(out) <- group
      class(out$n.studies.contributing) <-
        class(out$n.site) <-
        class(out$n.alt) <- "integer"

      ## 5b. Identify studies that actually contributed ------------------
      if (length(sv.list) > 0L) {
        effective_studies   <- which(!unlist(lapply(sv.list, is.null)))
        n_studies_effective <- length(effective_studies)
      } else {
        effective_studies   <- integer(0)
        n_studies_effective <- 0L
      }

      # Correct the number of contributing studies
      out$n.studies.contributing <- n_studies_effective

      ## 5c. Skip if no studies contributed ------------------------------
      if (n_studies_effective == 0L) {
        res <- dplyr::bind_rows(res, out)
        next
      }

      ## 5d. Sanity check: variant + covariance alignment per study ------
      for (i in effective_studies) {
        check1 <- sv.list[[i]][, "variant.id"] == colnames(V.list[[i]])
        check2 <- colnames(V.list[[i]])       == rownames(V.list[[i]])
        if (FALSE %in% c(check1, check2)) {
          stop(
            "Warning: for cohort ", i,
            " the variants in the single var file, colnames of covariance file,",
            " or rownames of covariance file do not match.\n"
          )
        }
      }

      ## 5e. Build global variant list for this transcript ---------------
      variant.list <- NULL
      for (i in effective_studies) {
        variant.list <- unique(c(variant.list, sv.list[[i]]$variant.id))
      }
      variant.list <- sort(variant.list)
      n.variants   <- length(variant.list)
      out$n.site   <- n.variants

      ## 5f. Reconstruct meta-analysis U and V ---------------------------
      U <- matrix(0, n.variants, 1,
                  dimnames = list(variant.list, "Score"))
      V <- matrix(0, n.variants, n.variants,
                  dimnames = list(variant.list, variant.list))

      for (i in effective_studies) {
        # Cohort-specific variant order (sorted)
        variant.list.cohort <- sv.list[[i]][, "variant.id"]
        variant.list.cohort <- sort(variant.list.cohort)

        # Ensure sv.list[[i]] rows are in the same order
        sv.list[[i]] <- sv.list[[i]][order(sv.list[[i]]$variant.id), ]

        # Add scores
        U[variant.list.cohort, ] <-
          U[variant.list.cohort, ] + sv.list[[i]][, "Score"]

        # Add full SPA-calibrated covariance (already aligned to cohort order)
        if (is.null(ncol(V.list[[i]]))) {
          V[variant.list.cohort, variant.list.cohort] <-
            matrix(V[variant.list.cohort, variant.list.cohort] + V.list[[i]])
        } else {
          V[variant.list.cohort, variant.list.cohort] <-
            matrix(
              V[variant.list.cohort, variant.list.cohort] +
                V.list[[i]][variant.list.cohort, variant.list.cohort]
            )
        }
      }

      ## 6. Run gene-based tests for this transcript ----------------------
      # Burden: used alone, and inside SMMAT
      U.sum <- sum(U[, "Score"])
      V.sum <- sum(V)
      GG1   <- rowSums(V)

      burden.pval <- pchisq(U.sum^2 / V.sum, df = 1, lower.tail = FALSE)
      out[, c("Burden.Score", "Burden.Variance", "Burden.pval")] <-
        c(U.sum, V.sum, burden.pval)
      class(out$Burden.Score) <-
        class(out$Burden.Variance) <-
        class(out$Burden.pval) <- "numeric"

      ## SKAT ------------------------------------------------------------
      if ("SKAT" %in% test) {
        Q <- sum(U^2)
        SKAT.pval        <- NA
        SKAT.pval.method <- NA
        if (mean(abs(V)) >= sqrt(.Machine$double.eps)) {
          pv <- regular(Q, V, n.variants)    # user-supplied function
          SKAT.pval        <- pv$pval
          SKAT.pval.method <- pv$method
        }
        out[, c("SKAT.pval", "SKAT.pval.method")] <- c(SKAT.pval, SKAT.pval.method)
        class(out$SKAT.pval)        <- "numeric"
        class(out$SKAT.pval.method) <- "character"
      }

      ## SKATO -----------------------------------------------------------
      if ("SKATO" %in% test) {
        Q <- sum(U^2)
        SKATO.pval        <- NA
        SKATO.pval.method <- NA
        if (mean(abs(V)) >= sqrt(.Machine$double.eps)) {
          res_skato <- GMMAT:::.skato_pval(U = U, V = V,
                                           rho = rho, method = "davies")
          Burden.Score    <- res_skato$Burden.score
          Burden.Variance <- res_skato$Burden.var
          Burden.pval     <- res_skato$Burden.pval
          SKAT.pval       <- res_skato$SKAT.pval
          SKATO.pval      <- res_skato$p
          SKATO.minp      <- res_skato$minp
          SKATO.minp.rho  <- res_skato$minp.rho
          out[, c("Burden.Score", "Burden.Variance", "Burden.pval",
                  "SKAT.pval", "SKATO.pval", "SKATO.minp", "SKATO.minp.rho")] <-
            c(Burden.Score, Burden.Variance, Burden.pval,
              SKAT.pval, SKATO.pval, SKATO.minp, SKATO.minp.rho)
        }
      }

      ## SMMAT -----------------------------------------------------------
      if ("SMMAT" %in% test) {
        # Burden-adjusted SKAT
        U2 <- U - GG1 * U.sum / V.sum
        Q  <- sum(U2^2)
        V2 <- V - tcrossprod(GG1) / V.sum

        theta.pval        <- NA
        theta.pval.method <- NA
        err               <- NA
        if (mean(abs(V2)) >= sqrt(.Machine$double.eps)) {
          pv <- regular(Q, V2, n.variants)
          theta.pval        <- pv$pval
          theta.pval.method <- pv$method
          err               <- pv$err
        }

        # Fisher’s method to combine burden and variance-component p-values
        SMMAT.pval <- tryCatch(
          pchisq(-2 * log(burden.pval) - 2 * log(theta.pval),
                 df = 4, lower.tail = FALSE),
          error = function(e) NA
        )

        if (is.na(SMMAT.pval)) {
          err        <- 1
          SMMAT.pval <- burden.pval
        }

        out[, c("theta.pval", "theta.pval.method", "err", "SMMAT.pval")] <-
          c(theta.pval, theta.pval.method, err, SMMAT.pval)
        class(out$theta.pval) <-
          class(out$err) <-
          class(out$SMMAT.pval) <- "numeric"
        class(out$theta.pval.method) <- "character"
      }

      ## 7. Append transcript-level result --------------------------------
      res <- dplyr::bind_rows(res, out)

    } # end transcript loop
  }   # end gene loop
      res<-data.frame(res,row.names=NULL)
  # close progress bar
  if (progress_bar && !is.null(pb)) close(pb)

  ## 8. Optional Cauchy meta p-value combination --------------------------
  if (combine.pval) {
    res <- comb_meta_pvalue(
      data    = res,
      pval.col = paste0(test, ".pval"),   # e.g. "Burden.pval", "SKAT.pval", ...
      min.cmac = min_meta_cmac
    )
  }

  return(res)
}



comb_meta_pvalue<-function(data,pval.col=c("Burden.pval","SKAT.pval","SKATO.pval","SMMAT.pval"),min.cmac=10){

result<-list()
pvalues<-names(data)[names(data) %in% pval.col]

data<-subset(data,n.alt>=min.cmac)

if (nrow(data) == 0) return(NULL)

if(length(pval.col)>1){
data$Cauchy.pval<-apply(data[,pvalues],1,function(x){x<-na.omit(x);CCT(x)})
}
genenames<-unique(data$Group)

sum1<-NULL
for (gg in 1:length(genenames)){
genename<-genenames[gg]
gdata<-subset(data,Group==genename & n.alt>=min.cmac)

if(nrow(gdata)>0){
gpvalue<-unlist(gdata[,pvalues])
gpvalue<-na.omit(gpvalue)
total.cauchy.pval<-CCT(gpvalue)
sum0<-data.frame(Group=genename,Cauchy.anytranscript.pval=total.cauchy.pval)
sum1<-rbind(sum1,sum0)
}
}

gres0<-subset(data,Transcript %in% sum1$Group)
gres1<-merge(gres0,sum1,by="Group")

result[["by.transcript"]]<-data
result[["combined.anytranscript"]]<-gres1

return(result)
}


###
###
CCT <- function(pvals, weights=NULL){
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }

  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }

  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  is.one <- (sum(pvals==1)>=1)
  if(is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if(is.zero){
    return(0)
  }
  if(is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }

  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }

  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }

  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}