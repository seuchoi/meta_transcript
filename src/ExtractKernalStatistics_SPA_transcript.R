####
#### 1) add SPA and 2)transcript analysis
testVariantSet_ExtractKernelStatistics_ScoresAndCovarianceMatrices_Sean <- function(nullmod, G, weights, var.info, neig = Inf, ntrace = Inf, Use.SPA=F, freq,
                                                                                       SAIGEGENEplus_collapse_threshold=1,grp){

       
                
        # Check for use.SPA, which is not yet supported
        #if(Use.SPA){
        #        stop("SPA not yet implemented for ExtractKernelStatistics function. Stopping.")
        #}

	# Modify var.info so output is in chr:pos:ref:alt format and can be compared across studies
        var.id.name <- paste0(var.info$chr, ":", var.info$pos, ":", var.info$ref, ":", var.info$alt)
        colnames(G) <- var.id.name


        if(is(G, "Matrix")){
                burden <- rowSums(G %*% Diagonal(x = weights))
                G <- G %*% Diagonal(x = weights)
        }else{
              	burden <- colSums(t(G) * weights)
                G <- t(t(G) * weights)
        }
        colnames(G) <- var.id.name # after tranformation colname disappeared! 02/01/2024


	if(is.null(nullmod$RSS0)){
                nullmod$RSS0 <- as.numeric(crossprod(nullmod$Ytilde))
        }

	# Calculate SKAT statistic
        U <- as.vector(crossprod(G, nullmod$resid)) # WGPY
        # SKAT test statistic
        Q <- sum(U^2)

        # adjust G for covariates and random effects
        burdentilde <- GENESIS:::calcGtilde(nullmod, burden)
        Gtilde <- GENESIS:::calcGtilde(nullmod, G) # P^{1/2}GW

        # Compute SKAT Variance
        ncolGtilde <- ncol(Gtilde)
        nrowGtilde <- nrow(Gtilde)

        if (ncolGtilde <= nrowGtilde) {
                V <- crossprod(Gtilde)
        } else {
                V <- tcrossprod(Gtilde)
        }

	    # We will start with single marker tests for each of the markers
        out <- GENESIS:::.testGenoSingleVarScore(Gtilde, G = G, resid = nullmod$resid, RSS0 = nullmod$RSS0)
        outnull <- out
        # Run SPA for each marker and the combined burden of very rare variants
        out <- SPA_pval_Sean(score.result = out, nullmod = nullmod, G = as.matrix(G), pval.thresh = 0.05)
        # Compute SPA adjusted variance
        out$SPA.Score.Variance <- (out$Score^2) / qchisq(out$SPA.pval, lower.tail=F, df=1)
        out[out$SPA.Score.Variance==0,'SPA.Score.Variance'] <- sqrt(outnull[which(out$SPA.Score.Variance==0),'Score.SE'])
        out[out$SPA.Score.Variance==0,'SPA.pval'] <- outnull[which(out$SPA.Score.Variance==0),'Score.pval']
        #out <- out[,c("Score", "SPA.Score.Variance", "SPA.pval", "Est", "Est.SE")]
        colnames(out)[c(5,6)] <- c("Raw.Est", "Raw.Est.SE")
        single_var_out <- out
        colnames(V)<-rownames(V)<-rownames(single_var_out)<-var.id.name # added

        # Compute SPA adjusted Sum of Variances (SAIGE-GENE, AJHG), we will use this for SKAT test later
        V_tilde <- single_var_out$SPA.Score.Variance
        Vsum_tilde <- sum(single_var_out$SPA.Score.Variance)
        
        # We will also compute a burden test for all markers
        out <- GENESIS:::.testGenoSingleVarScore(burdentilde, G = burden, resid = nullmod$resid, RSS0 = nullmod$RSS0)
        outnull <- out
        #Run SPA for the burden
        out <- SPA_pval_Sean(score.result = out, nullmod = nullmod, G = as.matrix(burden), pval.thresh = 0.05)
        # Compute SPA adjusted variance
        out$SPA.Score.Variance <- (out$Score^2) / qchisq(out$SPA.pval, lower.tail=F, df=1)
        out[out$SPA.Score.Variance==0,'SPA.Score.Variance'] <- sqrt(outnull[which(out$SPA.Score.Variance==0),'Score.SE'])
        out[out$SPA.Score.Variance==0,'SPA.pval'] <- outnull[which(out$SPA.Score.Variance==0),'Score.pval']
        #out <- out[,c("Score", "SPA.Score.Variance", "SPA.pval", "Est", "Est.SE")]
        colnames(out)[c(5,6)] <- c("Raw.Est", "Raw.Est.SE")
        colnames(out) <- paste0("Burden_", colnames(out))
        burden_out<-out
    
        # Compute SPA adjusted Burden Variance (SAIGE-GENE, AJHG)
        Vsum_downwardhat <- out$Burden_SPA.Score.Variance

        # Compute ratio to find more conservative variance (SAIGE-GENE, AJHG)
        r <- Vsum_tilde / Vsum_downwardhat
        r_tilde <- min(1, r)
        burden_out$r_tilde<-r_tilde

        # Adjust SKAT variance (SAIGE-GENE, AJHG)
        #diag(V) <- V_tilde
        #V <- V / r_tilde

        allvarlist<-grp[[1]]
        transcriptids<-unique(allvarlist$TranscriptID)
        av.transcriptids<-NULL
        
        
        ### run per transcript
        for (tp in 1:length(transcriptids)){
        transcriptid<-transcriptids[tp]
        subvarlist<-subset(allvarlist,TranscriptID==transcriptid)
        sub.var.id.name <- paste0(subvarlist@seqnames@values, ":", subvarlist@ranges@start, ":", subvarlist$ref, ":", subvarlist$alt)
        
        av.sub.var.id<-colnames(G)[colnames(G) %in% sub.var.id.name]
        if(length(av.sub.var.id)>0){
        av.transcriptids<-c(av.transcriptids,transcriptid)
        
        if(length(av.sub.var.id)==length(colnames(G))){
        
        new_burden_out<-burden_out[1,]
        burden_out<-rbind(burden_out,new_burden_out)    

        }else{
        
        newG<-G[,colnames(G) %in% av.sub.var.id]
        newweight<-weights[colnames(G) %in% av.sub.var.id]
        if(is(G, "Matrix")){
            newburden <- rowSums(newG %*% Diagonal(x = newweight))
            newG <- newG %*% Diagonal(x = newweight)
        }else{
          	newburden <- colSums(t(newG) * newweight)
            newG <- t(t(newG) * newweight)
        }

	    # Calculate SKAT statistic
        U <- as.vector(crossprod(newG, nullmod$resid)) # WGPY
        # SKAT test statistic
        Q <- sum(U^2)

        # adjust G for covariates and random effects
        burdentilde <- GENESIS:::calcGtilde(nullmod, newburden)
        Gtilde <- GENESIS:::calcGtilde(nullmod, newG) # P^{1/2}GW

        # Compute SKAT Variance
        ncolGtilde <- ncol(Gtilde)
        nrowGtilde <- nrow(Gtilde)

        if (ncolGtilde <= nrowGtilde) {
                newV <- crossprod(Gtilde)
        } else {
                newV <- tcrossprod(Gtilde)
        }
        
	    # We will start with single marker tests for each of the markers
        out <- GENESIS:::.testGenoSingleVarScore(Gtilde, G = newG, resid = nullmod$resid, RSS0 = nullmod$RSS0)
        outnull <- out
        # Run SPA for each marker and the combined burden of very rare variants
        out <- SPA_pval_Sean(score.result = out, nullmod = nullmod, G = as.matrix(newG), pval.thresh = 0.05)
        # Compute SPA adjusted variance
        out$SPA.Score.Variance <- (out$Score^2) / qchisq(out$SPA.pval, lower.tail=F, df=1)
        out[out$SPA.Score.Variance==0,'SPA.Score.Variance'] <- sqrt(outnull[which(out$SPA.Score.Variance==0),'Score.SE'])
        out[out$SPA.Score.Variance==0,'SPA.pval'] <- outnull[which(out$SPA.Score.Variance==0),'Score.pval']
        #out <- out[,c("Score", "SPA.Score.Variance", "SPA.pval", "Est", "Est.SE")]
        colnames(out)[c(5,6)] <- c("Raw.Est", "Raw.Est.SE")
        single_var_out <- out
        colnames(newV)<-rownames(newV)<-rownames(single_var_out)<-av.sub.var.id # added

        # Compute SPA adjusted Sum of Variances (SAIGE-GENE, AJHG), we will use this for SKAT test later
        newV_tilde <- single_var_out$SPA.Score.Variance
        newVsum_tilde <- sum(single_var_out$SPA.Score.Variance)
        
        # We will also compute a burden test for all markers
        out <- GENESIS:::.testGenoSingleVarScore(burdentilde, G = newburden, resid = nullmod$resid, RSS0 = nullmod$RSS0)
        outnull <- out
        #Run SPA for the burden
        out <- SPA_pval_Sean(score.result = out, nullmod = nullmod, G = as.matrix(newburden), pval.thresh = 0.05)
        # Compute SPA adjusted variance
        out$SPA.Score.Variance <- (out$Score^2) / qchisq(out$SPA.pval, lower.tail=F, df=1)
        out[out$SPA.Score.Variance==0,'SPA.Score.Variance'] <- sqrt(outnull[which(out$SPA.Score.Variance==0),'Score.SE'])
        out[out$SPA.Score.Variance==0,'SPA.pval'] <- outnull[which(out$SPA.Score.Variance==0),'Score.pval']
        #out <- out[,c("Score", "SPA.Score.Variance", "SPA.pval", "Est", "Est.SE")]
        colnames(out)[c(5,6)] <- c("Raw.Est", "Raw.Est.SE")
        colnames(out) <- paste0("Burden_", colnames(out))
        new_burden_out<-out
    
        # Compute SPA adjusted Burden Variance (SAIGE-GENE, AJHG)
        newVsum_downwardhat <- out$Burden_SPA.Score.Variance

        # Compute ratio to find more conservative variance (SAIGE-GENE, AJHG)
        newr <- newVsum_tilde / newVsum_downwardhat
        newr_tilde <- min(1, newr)
        new_burden_out$r_tilde<-newr_tilde
        burden_out<-rbind(burden_out,new_burden_out)
        }
        }else{}
        }
        burden_out$transcript<-c("all",av.transcriptids)
        burden_out$genename<-names(grp)
        out <- list(NULL)
        out[['burden_out']] <- burden_out
        out[['single_var_out']] <- single_var_out
        out[['covariance_matrix']] <- V
        return(out)
}


testVariantSet_Sean <- function( nullmod, G, weights, freq, use.weights=F, var.info,
                            test = c("Burden", "SKAT", "fastSKAT", "SMMAT", "fastSMMAT", "SKATO", "SKAT_SAIGEGENEplus", "ExtractKernelStatistics"),
                            burden.test = c("Score","Score.SPA"), collapse = FALSE, recessive=FALSE, recessive.model = c("strict", "putative"),
                            vc.type = "regular weighted", vc.test=c("Score","Score.SPA"), SAIGEGENEplus_collapse_threshold=10,grp=grp,
                            neig = 200, ntrace = 500,
                            rho = seq(from = 0, to = 1, by = 0.1)){
                           # pval.method = c("davies", "kuonen", "liu"),
                           # return.scores = FALSE, return.scores.cov = FALSE){

    test <- match.arg(test)
    burden.test <- match.arg(burden.test)
    vc.type <- match.arg(vc.type)
    # pval.method <- match.arg(pval.method)

    G <- GENESIS:::.genoAsMatrix(nullmod, G)
    if (test == "Burden") {
        if(collapse){
                burden.type <- "collapsing test"
        }else if(!use.weights){
                burden.type <- "regular burden"
        }else{
              	burden.type <- "externally weighted burden"
        }
	#cat('Running Burden test type', burden.type, 'using Pval method ', burden.test, '...\n')
        out <- testVariantSetBurden_Sean(nullmod, G, weights, burden.test = burden.test, collapse = collapse, recessive = recessive)
    }
    if (test == "SKAT") {
        if(vc.test=="Score.SPA"){
                stop('SPA is not yet implemented for', test, '...\n')
        }
	#cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSKAT_Sean(nullmod, G, weights, neig = Inf, ntrace = Inf)
                                   # return.scores, return.scores.cov)
    }
    if(test == "fastSKAT"){
        if(vc.test=="Score.SPA"){
                stop('SPA is not yet implemented for', test, '...\n')
        }
	#cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSKAT_Sean(nullmod, G, weights, neig, ntrace)
    }
    if (test == "SMMAT") {
        if(vc.test=="Score.SPA"){
                stop('SPA is not yet implemented for', test, '...\n')
        }
	#cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSMMAT_Sean(nullmod, G, weights, neig = Inf, ntrace = Inf)
    }
    if(test == "fastSMMAT"){
        if(vc.test=="Score.SPA"){
                stop('SPA is not yet implemented for', test, '...\n')
        }
	#cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSMMAT_Sean(nullmod, G, weights, neig, ntrace)
    }
    if(test == "SKATO"){
        if(vc.test=="Score.SPA"){
                stop('SPA is not yet implemented for', test, '...\n')
        }
	#cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSKATO_Sean(nullmod, G, weights, rho)
    }
    if(test == "SKAT_SAIGEGENEplus"){
        #cat('Running variance component-based test type', test, 'type', vc.type, 'using pvalue method', vc.test, '...\n')
        Use.SPA <- F
        if(vc.test=="Score.SPA"){
                Use.SPA <- T
        }
	out <- testVariantSetSKAT_SAIGEGENEplus_Sean(nullmod, G, weights, neig = Inf, ntrace = Inf, Use.SPA=Use.SPA, freq=freq,
                                                     SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold)
    }
    if(test == "ExtractKernelStatistics"){
        #cat('Extracting Kernel Statistics...\n')
        Use.SPA <- F
        if(vc.test=="Score.SPA"){
        	Use.SPA <- T
        }
	    out <- testVariantSet_ExtractKernelStatistics_ScoresAndCovarianceMatrices_Sean(nullmod, G, weights, var.info, neig = Inf, ntrace = Inf, Use.SPA=Use.SPA, freq=freq,
                                                                                       SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold,grp=grp)
    }
    return(out)
}






setGeneric("assocTestAggregate_Sean", function(gdsobj, ...) standardGeneric("assocTestAggregate_Sean"))

match.arg_Sean <- function(test) {
    if (length(test) > 1) test <- NULL
    match.arg(test, choices=c("Burden", "SKAT", "fastSKAT", "SMMAT", "fastSMMAT", "SKATO", "SKAT_SAIGEGENEplus", "ExtractKernelStatistics"))
}

setMethod("assocTestAggregate_Sean",
          "SeqVarIterator",
          function(gdsobj, null.model, AF.max=1, MAC.max=Inf, use.weights=F,
                   weight.beta=c(1,1), weight.user=NULL,
                   test=c("Burden", "SKAT", "fastSKAT", "SMMAT", "SKATO", "SKAT_SAIGEGENEplus", "ExtractKernelStatistics"),
                   burden.test=c("Score", "Score.SPA"), collapse=FALSE, recessive=F, recessive.model=c("strict", "putative"),
                   vc.test=c("Score", "Score.SPA"), vc.type="regular weighted", SAIGEGENEplus_collapse_threshold=10,grp=gr,
                   # pval.method=c("davies", "kuonen", "liu"),
                   neig = 200, ntrace = 500,
                   rho = seq(from = 0, to = 1, by = 0.1),
                   sparse=TRUE, imputed=FALSE,
                   male.diploid=TRUE, genome.build=c("hg19", "hg38"),
                   verbose=TRUE) {

              # check argument values
              test <- match.arg_Sean(test)
              burden.test <- match.arg(burden.test)
              # pval.method <- match.arg(pval.method)

              # don't use sparse matrices for imputed dosages
              if (imputed) sparse <- FALSE

              # coerce null.model if necessary
              if (sparse) null.model <- GENESIS:::.nullModelAsMatrix(null.model)

              # filter samples to match null model
              sample.index <- GENESIS:::.setFilterNullModel(gdsobj, null.model, verbose=verbose)

              # do we need to match on alleles?
              match.alleles <- any(c("ref", "alt") %in% names(mcols(currentRanges(gdsobj))))

              # check ploidy
              if (SeqVarTools:::.ploidy(gdsobj) == 1) male.diploid <- FALSE

              # results
              res <- list()
              res.var <- list()
              if(test == "ExtractKernelStatistics"){
                  res.covariance <- list()
              }

              i <- 1
              n.iter <- length(variantFilter(gdsobj))
              set.messages <- ceiling(n.iter / 100) # max messages = 100
              iterate <- TRUE
              while (iterate) {
                  var.info <- variantInfo(gdsobj, alleles=match.alleles, expanded=TRUE)

                  if (!imputed) {
                      geno <- expandedAltDosage(gdsobj, use.names=FALSE, sparse=sparse)[sample.index,,drop=FALSE]
                  } else {
                      geno <- imputedDosage(gdsobj, use.names=FALSE)[sample.index,,drop=FALSE]
                  }

                  if (match.alleles) {
                      index <- GENESIS:::.matchAlleles(gdsobj, var.info)
                      var.info <- var.info[index,,drop=FALSE]
                      geno <- geno[,index,drop=FALSE]
                  } else {
                      index <- NULL
                  }

                  # number of non-missing samples
                  # n.obs <- colSums(!is.na(geno))
                  n.obs <- GENESIS:::.countNonMissing(geno, MARGIN = 2)

                  # allele frequency
                  freq <- GENESIS:::.alleleFreq(gdsobj, geno, variant.index=index, sample.index=sample.index,
                                      male.diploid=male.diploid, genome.build=genome.build)
                  #freq <- GENESIS:::.alleleFreq(geno) ## added Oct/6/2023
                  # filter monomorphic variants
                  keep <- GENESIS:::.filterMonomorphic(geno, count=n.obs, freq=freq$freq, imputed=imputed)

                  # exclude variants with freq > max & MAC > max
                  keep <-  keep & freq$freq <= AF.max & freq$MAC <= MAC.max
                  if (!all(keep)) {
                      var.info <- var.info[keep,,drop=FALSE]
                      geno <- geno[,keep,drop=FALSE]
                      n.obs <- n.obs[keep]
                      freq <- freq[keep,,drop=FALSE]
                  }

                  # weights
                  if (is.null(weight.user)) {
                      # Beta weights
                      weight <- GENESIS:::.weightFromFreq(freq$freq, weight.beta)
                  } else {
                      # user supplied weights
                      weight <- currentVariants(gdsobj)[[weight.user]][expandedVariantIndex(gdsobj)]
                      if (!is.null(index)) weight <- weight[index]
                      weight <- weight[keep]

                      weight0 <- is.na(weight) | weight == 0
                      if (any(weight0)) {
                          keep <- !weight0
                          var.info <- var.info[keep,,drop=FALSE]
                          geno <- geno[,keep,drop=FALSE]
                          n.obs <- n.obs[keep]
                          freq <- freq[keep,,drop=FALSE]
                          weight <- weight[keep]
                      }
                  }

                  # number of variant sites
                  n.site <- length(unique(var.info$variant.id))

                  # number of alternate alleles
                  n.alt <- sum(geno, na.rm=TRUE)

                  # number of samples with observed alternate alleles > 0
                  n.sample.alt <- sum(rowSums(geno, na.rm=TRUE) >= 0.5)

                  # creat variant ids
                  var.info$var.id<-paste(var.info$chr,var.info$pos,var.info$ref,var.info$alt,sep=":")
                  
                  ## trnscript
                  allvarlist<-grp[[i]]
                  transcriptids<-unique(allvarlist$TranscriptID)
                  av.transcriptids<-NULL        
                  
                  ### run per transcript
                    for (tp in 1:length(transcriptids)){
                    #print(tp)
                    transcriptid<-transcriptids[tp]
                    subvarlist<-subset(allvarlist,TranscriptID==transcriptid)
                    sub.var.id.name <- paste0(subvarlist@seqnames@values, ":", subvarlist@ranges@start, ":", subvarlist$ref, ":", subvarlist$alt)
                    av.sub.var.id<- var.info$var.id[var.info$var.id %in% sub.var.id.name]
                    colnums<-which(var.info$var.id %in% sub.var.id.name)
                    # number of variant sites
                    n.site1<-length(av.sub.var.id)
                    n.site<-c(n.site,n.site1)
                    # number of alternate alleles
                    n.alt1<- sum(geno[,colnums,drop=FALSE], na.rm=TRUE)
                    n.alt<-c(n.alt,n.alt1)
                    # number of samples with observed alternate alleles > 0
                    n.sample.alt1 <- sum(rowSums(geno[,colnums,drop=FALSE], na.rm=TRUE) >= 0.5)
                    n.sample.alt<-c(n.sample.alt,n.sample.alt1)
                    }
                    
                    # keep the site>0 transcript
                    countout<-data.frame(n.site, n.alt, n.sample.alt)
                    avtranscripts<-transcriptids[which(n.site[2:length(n.site)]>0)]    
                    grp[[i]]<-subset(grp[[i]],TranscriptID %in% avtranscripts)
                    res[[i]] <-subset(countout,n.site>0) 
                    res.var[[i]] <- cbind(var.info, n.obs, freq, weight)
                  if(test == "ExtractKernelStatistics"){
                      cat('ExtractKernelStatistics number', i, '\n')
                      res.covariance[[i]] <- NA
                  }
		  not_run <- FALSE
                  if (n.site[1] > 0) {
                      # mean impute missing values, unless it is collapsing test in which case we will impute to zero
		      if(collapse){
                          if (any(n.obs < nrow(geno))) {
                                geno <- zeroImpute_Sean(geno, freq$freq)
                          }
                      }else{
                          if (any(n.obs < nrow(geno))) {
                                geno <- meanImpute_Sean(geno, freq$freq)
                          }
                      }

		      # if strict recessive analysis code for that
		      if(recessive){
			  if(recessive.model=="strict"){
		          	geno <- recessive_strict_coding_Sean(geno)
			  	n.alt <- sum(geno, na.rm=TRUE)
			  	n.sample.alt <- sum(rowSums(geno, na.rm=TRUE) >= 0.75)
			  }else{
				geno <- geno/2
			        n.alt <- n.sample.alt <- sum(rowSums(geno, na.rm=TRUE) >= 0.75)
			  }
			  res[[i]][2] <- n.alt
			  res[[i]][3] <- n.sample.alt
			  if(n.alt==0){
				  not_run <- TRUE
			  }
		      }

                      if(!not_run){
			   #do the test
			   assoc <- testVariantSet_Sean(null.model, G=geno, use.weights=use.weights, weights=weight, freq=freq,
                                              test=test, burden.test=burden.test, collapse=collapse, recessive=recessive, recessive.model=recessive.model,
					      var.info=var.info,
                                              vc.test=vc.test, vc.type=vc.type, SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold,
                                              neig = neig, ntrace = ntrace,
                                              rho=rho,grp=grp[i])
                                              # pval.method=pval.method)
                      	   if(test == 'ExtractKernelStatistics'){
                           	     	res[[i]] <- cbind(res[[i]], assoc[['burden_out']], stringsAsFactors=FALSE)
                                	res.var[[i]]$variant.id <- paste0(res.var[[i]]$chr, ":", res.var[[i]]$pos, ":", res.var[[i]]$ref, ":", res.var[[i]]$alt)
                                	assoc[['single_var_out']]$variant.id <- rownames(assoc[['single_var_out']])
                                	res.var[[i]] <- merge(res.var[[i]], assoc[['single_var_out']], by="variant.id", all=T)
                                	res.covariance[[i]] <- assoc[['covariance_matrix']]
                      	   }else{
                            		res[[i]] <- cbind(res[[i]], assoc, stringsAsFactors=FALSE)
                      	   }
		      }
                  }

                  if (verbose & n.iter > 1 & i %% set.messages == 0) {
                      message(paste("Iteration", i , "of", n.iter, "completed"))
                  }
                  i <- i + 1
                  iterate <- SeqVarTools:::iterateFilter(gdsobj, verbose=F)
              }
              if(test == 'ExtractKernelStatistics'){
                  res <- list(results=dplyr::bind_rows(res), variantInfo=res.var, covariance_matrix=res.covariance)
                  names(res$variantInfo) <- names(grp)
                  names(res$covariance_matrix) <- names(grp)
                  out_res<-res
              }else{
                  res <- list(results=dplyr::bind_rows(res), variantInfo=res.var)
                  out_res <- GENESIS:::.annotateAssoc(gdsobj, res)
              }
              return(out_res)
          })



kernell_variance_component_v2<-function(gdsfile, groupfile, phenfile, ID_col, nullfile, outfile,
                                       AF.max=0.001, MAC.max=Inf, use.weights=FALSE,
                                       vc.test=c("Score", "Score.SPA"),
                                       test=c("SKAT", "SKATO", "SMMAT", "SKAT_SAIGEGENEplus", "ExtractKernelStatistics"),
                                       SAIGEGENEplus_collapse_threshold=10, weight.beta=c(1,1)){
        #'
        #' gdsfile = string specifying the file name of the genetic dataset; dataset should be in SeqArray GDS format
        #' groupfile = string specifyinng the file name of the grouping file; the grouping file contains information of variants to be included in the analysis:
        #'             The grouping file should be a single dataframe called 'group' that is saved within a .RData file
        #'             The dataframe should contain the following columns in this order: varid, group_id, chr, pos, ref, alt. All other columns are optional.
        #'             Optionally, a column named 'weight' can be added for weighted burden tests.
        #'             An example of a grouping dataframe bellow:
        #'
        #'                           varid        group_id chr       pos ref alt         func Dscore
        #'             1 1:100007074:CTG:C ENSG00000283761   1 100007074 CTG   C hclof_noflag     NA
        #'             2 1:100007074:CTG:C ENSG00000117620   1 100007074 CTG   C hclof_noflag     NA
        #'             3   1:100007098:T:C ENSG00000283761   1 100007098   T   C     missense     26
        #'             4   1:100007098:T:C ENSG00000117620   1 100007098   T   C     missense     26
        #'             5   1:100007109:C:T ENSG00000283761   1 100007109   C   T hclof_noflag     NA
        #'             6   1:100007109:C:T ENSG00000117620   1 100007109   C   T hclof_noflag     NA
        #'               Dtools    Weight gnomAD_AFR_AMR_EAS_NFE_SAS_POPMAX
        #'             1     NA 1.0000000                                 0
        #'             2     NA 1.0000000                                 0
        #'             3     28 0.9285714                                 0
        #'             4     28 0.9285714                                 0
        #'             5     NA 1.0000000                                 0
        #'             6     NA 1.0000000                                 0
        #'
        #' phenfile = string specifying the phenotype file; phenotype file should be in .tsv format.
        #'            Phenotype file should contain sample identifiers (that match those in the GDS file), the outcome variable, and any fixed-effects covariates.
        #' ID_col = string specifying the column name for the column containing the sample ID information
        #' nullfile = string specifying the null-model file; this file contains the null-model that can be made using the 'fitNullModel' function from GENESIS or using our fit_nullmodel function.
        #' outfile = string specifying the preferred output location for the gene-based results.
        #' AF.max = numeric specifying the maximum allele frequency for including variants in the analysis. Variants with MAF>AF.max will be removed.
        #' MAC.max = numeric specifying the maximum minor allele count for including variants in the analysis. Variants with MAC>MAC.max will be removed.
        #' use.weights = logical indicating whether to use external weights in the burden test. Only works for collapse = FALSE. A column called 'weight' should be included in the grouping file.
        #' vc.test = vector of kernell-based tests to perform.


        if("Burden" %in% test){
                stop("Burden type test is not supported by this function. For burden use 'hclofburden()'. Stopping run.")
        }

        if(use.weights==F){
                vc.type <- "regular weighted"
        }else{
                vc.type <- "externally weighted"
                weight.beta <- c(1,1)
                cat("Note: because weights are pre-specified, the c(1,1) beta distribution (uniform distribution) will be used.\n")
        }

        cat(paste0('\n\nVariance component test type is ', vc.type, ' ', test, ' using pvalue method ', vc.test, ' with beta distribution of ', paste0("(", weight.beta[1], ",", weight.beta[2], ")"), '.\n\n\n'))

        # Samples
        phen1<-fread(phenfile,header=T,data.table=F,sep="\t")
        names(phen1)[which(colnames(phen1)==ID_col)]<-"sample.id"
        id_int <- FALSE
        if(class(phen1$sample.id)=='integer'){
                id_int <- TRUE
                class(phen1$sample.id) <- 'character'
        }
        samid0<-phen1$sample.id

        # Read gds file
        gds <- seqOpen(gdsfile, allow.duplicate=T)
        samples <- seqGetData(gds, "sample.id")
        if(id_int){class(samples)<-"character"}
        missamples<-samples[!samples %in% samid0]
        misphen<-data.frame(matrix(NA,nrow=length(missamples),ncol=ncol(phen1)))
        colnames(misphen)<-names(phen1)
        misphen$sample.id<-missamples
        combphen<-rbind(phen1,misphen)
        rownames(combphen)<-combphen$sample.id
        combphen2<-combphen[samples,]
        #if(id_int){class(combphen2$sample.id) <- 'integer'}

        # Construct a SeqVarData object
        seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(combphen2))

        # Filter the gdsfile
        seqSetFilter(seqData, sample.id=samid0)

        # Annotation file
        annot<-get(load(groupfile))
        annot <- as.data.frame(annot)
        #class(annot$chr) <- "numeric"
        class(annot$pos) <- "numeric"

        # Grouping file; add weights if weights are selected
        weights.found<-FALSE
        if(use.weights){
                if(!"weight" %in% colnames(annot)){
                        cat("\nWARNING: no column named 'weight' found in the grouping file; no weights will be applied.\n")
                        gr<-aggregateGRangesList(annot)
                }else{
                        #annot <- annot[,c("group_id", "chr", "pos", "ref", "alt", "weight")]
                        cat("\nuse.weights=T and 'weight' column found in grouping file; variant weights will be applied.\n")
                        gr<-aggregateGRangesList(annot)
                        weights.found<-TRUE
                }
        }else{
                gr<-aggregateGRangesList(annot)
        }

        # Create the iterator
        iterator <- SeqVarListIterator(seqData, variantRanges=gr)

        # Load null model
        nullmod<-get(load(nullfile))

        # Perfrom assocation test; apply weights if provided
        if(weights.found){
                assoc <- assocTestAggregate_Sean(iterator, nullmod, AF.max=AF.max, MAC.max=MAC.max, test=test, vc.test=vc.test, vc.type=vc.type, collapse = FALSE, verbose=TRUE, use.weights=T, weight.user="weight",
                                                 SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold, weight.beta=c(1,1),gr=gr)
        }else{
                assoc <- assocTestAggregate_Sean(iterator, nullmod, AF.max=AF.max, MAC.max=MAC.max, test=test, vc.test=vc.test, vc.type=vc.type, collapse = FALSE, verbose=TRUE, use.weight=F,
                                                 SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold, weight.beta=weight.beta,gr=gr)
        }

        # Save results
        save(assoc,file=outfile)
        seqClose(gds)
}



kernell_variance_component_aou<-function(gdsfile, groupfile, phenfile, ID_col, nullfile, outfile,
                                       AF.max=0.001, MAC.max=Inf, use.weights=FALSE,
                                       vc.test=c("Score", "Score.SPA"),
                                       test=c("SKAT", "SKATO", "SMMAT", "SKAT_SAIGEGENEplus", "ExtractKernelStatistics"),
                                       SAIGEGENEplus_collapse_threshold=10, weight.beta=c(1,1)){
        #'
        #' gdsfile = string specifying the file name of the genetic dataset; dataset should be in SeqArray GDS format
        #' groupfile = string specifyinng the file name of the grouping file; the grouping file contains information of variants to be included in the analysis:
        #'             The grouping file should be a single dataframe called 'group' that is saved within a .RData file
        #'             The dataframe should contain the following columns in this order: varid, group_id, chr, pos, ref, alt. All other columns are optional.
        #'             Optionally, a column named 'weight' can be added for weighted burden tests.
        #'             An example of a grouping dataframe bellow:
        #'
        #'                           varid        group_id chr       pos ref alt         func Dscore
        #'             1 1:100007074:CTG:C ENSG00000283761   1 100007074 CTG   C hclof_noflag     NA
        #'             2 1:100007074:CTG:C ENSG00000117620   1 100007074 CTG   C hclof_noflag     NA
        #'             3   1:100007098:T:C ENSG00000283761   1 100007098   T   C     missense     26
        #'             4   1:100007098:T:C ENSG00000117620   1 100007098   T   C     missense     26
        #'             5   1:100007109:C:T ENSG00000283761   1 100007109   C   T hclof_noflag     NA
        #'             6   1:100007109:C:T ENSG00000117620   1 100007109   C   T hclof_noflag     NA
        #'               Dtools    Weight gnomAD_AFR_AMR_EAS_NFE_SAS_POPMAX
        #'             1     NA 1.0000000                                 0
        #'             2     NA 1.0000000                                 0
        #'             3     28 0.9285714                                 0
        #'             4     28 0.9285714                                 0
        #'             5     NA 1.0000000                                 0
        #'             6     NA 1.0000000                                 0
        #'
        #' phenfile = string specifying the phenotype file; phenotype file should be in .tsv format.
        #'            Phenotype file should contain sample identifiers (that match those in the GDS file), the outcome variable, and any fixed-effects covariates.
        #' ID_col = string specifying the column name for the column containing the sample ID information
        #' nullfile = string specifying the null-model file; this file contains the null-model that can be made using the 'fitNullModel' function from GENESIS or using our fit_nullmodel function.
        #' outfile = string specifying the preferred output location for the gene-based results.
        #' AF.max = numeric specifying the maximum allele frequency for including variants in the analysis. Variants with MAF>AF.max will be removed.
        #' MAC.max = numeric specifying the maximum minor allele count for including variants in the analysis. Variants with MAC>MAC.max will be removed.
        #' use.weights = logical indicating whether to use external weights in the burden test. Only works for collapse = FALSE. A column called 'weight' should be included in the grouping file.
        #' vc.test = vector of kernell-based tests to perform.


        if("Burden" %in% test){
                stop("Burden type test is not supported by this function. For burden use 'hclofburden()'. Stopping run.")
        }

        if(use.weights==F){
                vc.type <- "regular weighted"
        }else{
                vc.type <- "externally weighted"
                weight.beta <- c(1,1)
                cat("Note: because weights are pre-specified, the c(1,1) beta distribution (uniform distribution) will be used.\n")
        }

        cat(paste0('\n\nVariance component test type is ', vc.type, ' ', test, ' using pvalue method ', vc.test, ' with beta distribution of ', paste0("(", weight.beta[1], ",", weight.beta[2], ")"), '.\n\n\n'))

        # Samples
        phen1<-fread(phenfile,header=T,data.table=F,sep="\t")
        names(phen1)[which(colnames(phen1)==ID_col)]<-"sample.id"
        id_int <- FALSE
        if(class(phen1$sample.id)=='integer'){
                id_int <- TRUE
                class(phen1$sample.id) <- 'character'
        }
        samid0<-phen1$sample.id

        # Read gds file
        gds <- seqOpen(gdsfile, allow.duplicate=T)
        samples <- seqGetData(gds, "sample.id")
        if(id_int){class(samples)<-"character"}
        missamples<-samples[!samples %in% samid0]
        misphen<-data.frame(matrix(NA,nrow=length(missamples),ncol=ncol(phen1)))
        colnames(misphen)<-names(phen1)
        misphen$sample.id<-missamples
        combphen<-rbind(phen1,misphen)
        rownames(combphen)<-combphen$sample.id
        combphen2<-combphen[samples,]
        if(id_int){class(combphen2$sample.id) <- 'integer'}

        # Construct a SeqVarData object
        seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(combphen2))

        # Filter the gdsfile
        seqSetFilter(seqData, sample.id=samid0)

        # Annotation file
        annot<-get(load(groupfile))
        annot <- as.data.frame(annot)[1:10]
        #class(annot$chr) <- "numeric"
        class(annot$pos) <- "numeric"

        # Grouping file; add weights if weights are selected
        weights.found<-FALSE
        if(use.weights){
                if(!"weight" %in% colnames(annot)){
                        cat("\nWARNING: no column named 'weight' found in the grouping file; no weights will be applied.\n")
                        gr<-aggregateGRangesList(annot)
                }else{
                        #annot <- annot[,c("group_id", "chr", "pos", "ref", "alt", "weight")]
                        cat("\nuse.weights=T and 'weight' column found in grouping file; variant weights will be applied.\n")
                        gr<-aggregateGRangesList(annot)
                        weights.found<-TRUE
                }
        }else{
                gr<-aggregateGRangesList(annot)
        }

        # Create the iterator
        iterator <- SeqVarListIterator(seqData, variantRanges=gr)

        # Load null model
        nullmod<-get(load(nullfile))

        # Perfrom assocation test; apply weights if provided
        if(weights.found){
                assoc <- assocTestAggregate_Sean(iterator, nullmod, AF.max=AF.max, MAC.max=MAC.max, test=test, vc.test=vc.test, vc.type=vc.type, collapse = FALSE, verbose=TRUE, use.weights=T, weight.user="weight",
                                                 SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold, weight.beta=c(1,1),gr=gr)
        }else{
                assoc <- assocTestAggregate_Sean(iterator, nullmod, AF.max=AF.max, MAC.max=MAC.max, test=test, vc.test=vc.test, vc.type=vc.type, collapse = FALSE, verbose=TRUE, use.weight=F,
                                                 SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold, weight.beta=weight.beta,gr=gr)
        }

        # Save results
        save(assoc,file=outfile)
        seqClose(gds)
}
