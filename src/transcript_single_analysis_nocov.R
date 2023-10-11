###
### meta-analysis function
### Function 1

aa0<-transcript_single_analysis_nocovout(study_path="/home/jupyter/workspaces/cadrvas/result/cad_185k_CAD_chr22.RData",
grouping_path="/home/jupyter/workspaces/cadrvas/annotation/AoU_250K_exome_annotation_VEP105_chr22.vcf.gz.hclof_noflag_POPMAX0.001.RData",test="Burden",
min_study_cmac=1,min_meta_cmac=20,use.anytranscript=T,combine.pval=T)

transcript_single_analysis_nocovout<-function(study_path=study_path,
                                   grouping_path=grouping_path,test=c("Burden"),
                                   min_study_cmac=1,min_meta_cmac=20,
                                   use.anytranscript=T,combine.pval=T){ #(0)#

### output filename
### grouping filename: transcript_id column is required
### any transcript is optional
### available tests

## read and combining results from 22 chromosomes for each phenotype
   n_studies  <- length(study_path)
   study_list <- list()
   for(i in 1:n_studies){
       load(study_path[i])
       study_list[[i]] <- assoc
   }

## Find number of unique SNP groups (by GENE) from results
## ver 2
   groupings <- lapply(study_list, function(x) rownames(x$results) )
   groupings <- do.call(c, groupings)
   groupings <- unique(groupings)
   length(groupings)   # 1081

## Determine grouping label group_id=(ENSG...)  transcript_id=(ENST...)
   transgroup_list<-list()
   for(i in 1:n_studies){
     if(use.anytranscript){ # then group_id = c(ENSG..., ENST...)
        grp1<-grp0<-get(load(grouping_path[[i]]))
        names(grp1)[which(names(grp1) %in% c("TranscriptID"))]<-names(grp0)[which(names(grp0) %in% c("TranscriptID"))]<-"transcript_id"
        grp1$transcript_id  <-grp1$group_id
        transgroup_list[[i]]<-rbind(grp0,grp1)
        grp0 <- NULL
        grp1 <- NULL
     }else{
        transgroup_list[[i]]<-get(load(grouping_path[[i]]))
     }
   }

## ---------------------------------------------------------------------
## Cycle over Gene groupings
   res <- NULL
   num <- 1
   n__groupings <- length(groupings)


## start gene level
   for(group in groupings){#(1)#
    outnums <- c(1,seq(0,n__groupings,by=10),n__groupings)    # ???????
    print(num)
    if(num %in% outnums){
       cat('Busy with group', group, 'which is', num, 'out of', n__groupings, '...\n')
    }
    num <- num + 1

    # Find transcript group
    transcripts <- lapply(transgroup_list, function(x) unique(x$transcript_id[which(x$group_id == group)]) )
    transcripts <- unique(do.call(c, transcripts))
    transcripts <- transcripts[!is.na(transcripts)]

    if(length(transcripts)>0){ #(2)#
      for(tt in 1:length(transcripts)){ #(3)#
         print(tt)
         transID <- transcripts[tt]
         n_studies_effective <- NA
         variant.id <- n.alts <- NULL
         V.list <- sv.list <- list()

      # --------------------------------------------------------------------
      # creating sv.list and V.list
      for(i in 1:n_studies){   #(4)#

            transgroup_data      <- transgroup_list[[i]]
            trans_group_variants <- subset(transgroup_data,transcript_id==transID)

          if(nrow(trans_group_variants)>0){
          # Generate ChrPosRefAlt from Annotation
            trans_group_variants$variant.id<-paste(trans_group_variants$chr,trans_group_variants$pos,
                                                   trans_group_variants$ref,trans_group_variants$alt,sep=":")

          # info from model results
            total_variant <- study_list[[i]]$variantInfo[[group]]
            total_cov     <- study_list[[i]]$covariance_matrix[[group]]

          # overlap between results and annotation
            sub_score <- total_variant[total_variant$variant.id %in% trans_group_variants$variant.id,]
            sub_score <- unique(sub_score)

          # update variant.id, n.alts, sv.list, V.list only if sum(MAC) >= threshold
            if(sum(sub_score$MAC)>=min_study_cmac){
              variant.id <- c(variant.id, sub_score$variant.id)
              n.alts     <- c(n.alts, sum(sub_score$MAC))

              sub_score    <- sub_score[order(sub_score$variant.id),]
              sv.list[[i]] <- sub_score
              vindex       <- sub_score$variant.id
              V.list[[i]]  <- total_cov[vindex,vindex]
            }
          }
          # creating summary file "out"  ??????????????????????????????????????????????????
          out <- data.frame(group, transcript=transID, n_studies_effective=length(sv.list),
                            n.site=length(unique(variant.id)), n.alt=sum(n.alts), stringsAsFactors=F)
          colnames(out) <- c("Group", "Transcript","n.studies.contributing", "n.site", "n.alt")
          rownames(out) <- group
          class(out$n.studies.contributing) <- class(out$n.site) <- class(out$n.alt) <- "integer"
        } #(4)#

      # Fix when one sv list is 0
        if(length(sv.list)>0){
           effective_studies   <- which(!unlist(lapply(sv.list,is.null)))
           n_studies_effective <- length(effective_studies)
        }else{
           effective_studies   <- NA
           n_studies_effective <- 0
        }
      # --------------------------------------------------------------------

      # --------------------------------------------------------------------
        if(n_studies_effective>0){ #(5)#
         # check
           for(i in c(effective_studies)){
              check  <- nrow(sv.list[[i]])
              check1 <- sv.list[[i]]$variant.id == colnames(V.list[[i]])
              check2 <- colnames(V.list[[i]]) == rownames(V.list[[i]])
              if(F %in% c(check1, check2)){
                stop("Warning: for cohort ", i,
                     "the variants in the single var file, colnames of Covariance file, or rownames of covariance file do not match.\n")
              }
           }

         # Prepare
           variant.list <- lapply(sv.list, function(x) x$variant.id)
           variant.list <- unique(do.call(c, variant.list))
           variant.list <- variant.list[order(variant.list)]

           n.site <- n.variants <- length(variant.list)
           out$n.site <- n.site
           print(length(variant.list))

         # Reconstruct variant and covariance matrix
           U <- matrix(0, n.variants, 1, dimnames=list(variant.list, 'Score'))
           V <- matrix(0, n.variants, n.variants, dimnames=list(variant.list, variant.list))

           for(i in c(effective_studies)){
             variant.list.cohort     <- sv.list[[i]]$variant.id
             variant.list.cohort     <- variant.list.cohort[order(variant.list.cohort)]
             sv.list[[i]]            <- sv.list[[i]][order(sv.list[[i]]$variant.id),]
             U[variant.list.cohort,] <- U[variant.list.cohort,] + sv.list[[i]]$Score
             if(is.null(ncol(V.list[[i]]))){
               V[variant.list.cohort, variant.list.cohort] <- matrix(V[variant.list.cohort, variant.list.cohort] + V.list[[i]])
             }else{
               V[variant.list.cohort, variant.list.cohort] <- matrix(V[variant.list.cohort, variant.list.cohort] + V.list[[i]][variant.list.cohort,variant.list.cohort])
             }
           }

         #############################
         ###        Run test       ###
         #############################
         # Start with Burden, this is also used for SMMAT test, and also for SKAT when implemented later
           U.sum <- sum(U[,'Score'])
           V.sum <- sum(V)
#           GG1   <- rowSums(V)

         # Then using adapted GENESIS script
         ### Burden test

            if("Burden" %in% test){
              out0<-burden.fun(U.sum=U.sum,V.sum=V.sum)
             out[,names(out0)]<-unlist(out0)
           }
         #############################

        } #(5)#
        res <- dplyr::bind_rows(res, out)

      } #(3)#
    } #(2)#

   } #(1)#

## combine p-values using cauchy distribution
   if(combine.pval){
     res<-comb_meta_pvalue(data=res,min.cmac=min_meta_cmac)
   }

return(res)
} #(0)#


####
#### combined p-value from meta-analysis

comb_meta_pvalue<-function(data,pval.col=c("Burden.pval","SKAT.pval","SKATO.pval","SMMAT.pval"),min.cmac=10){

result<-list()
pvalues<-names(data)[names(data) %in% pval.col]

data<-subset(data,n.alt>=min.cmac)
if(length(pvalues)!=1){
  data$Cauchy.pval<-apply(data[,pvalues],1,function(x){x<-na.omit(x);CCT(x)})
}else{
  data$Cauchy.pval<-data[,pvalues]
}
genenames<-unique(data$Group)

sum1<-NULL
for (gg in 1:length(genenames)){
genename<-genenames[gg]
gdata<-subset(data,Group==genename &n.alt>=min.cmac)

if(nrow(gdata)>0){
gpvalue<-unique(unlist(gdata[,pvalues]))
gpvalue<-na.omit(gpvalue)
if(length(gpvalue)!=1){
total.cauchy.pval<-CCT(gpvalue)
}else{
total.cauchy.pval<-gpvalue
}

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
