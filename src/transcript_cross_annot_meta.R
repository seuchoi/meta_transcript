###
###
### perfrom meta-analysis of
#study_path<-c("test1.RData","test2.RData")
                             
transcript_cross_annot_meta<-function(study_path=study_path,test=c("Burden"),
                                   min_study_cmac=1,min_meta_cmac=20,
                                   use.anytranscript=T){ #(0)#

### output filename
### grouping filename: transcript_id column is required
### any transcript is optional
### available tests

## read and combining results from 22 chromosomes for each phenotype
   n_studies  <- length(study_path)
   study_list <- list()
   for(i in 1:n_studies){
       assoc<-get(load(study_path[i]))
       study_list[[i]] <- assoc
   }

## Find number of unique SNP groups (by GENE) from results
## ver 2
   groupings <- lapply(study_list, function(x) unique(x$Group) )
   groupings <- do.call(c, groupings)
   groupings <- unique(groupings)
   length(groupings)   # 1081

## ---------------------------------------------------------------------
## Cycle over Gene groupings
   res <- NULL
   num <- 1
   n__groupings <- length(groupings)

   result<-list()
   cauchy.result<-NULL
## start gene level
##   for(group in groupings[1:10]){#(1)#
    for (gg in 1:length(groupings)){
    group<-groupings[gg]
    outnums <- c(1,seq(0,n__groupings,by=10),n__groupings)    # ???????
    #print(num)
    if(num %in% outnums){
       cat('Busy with group', group, 'which is', num, 'out of', n__groupings, '...\n')
    }
    num <- num + 1

    gpvalue<-NULL
    pvalues<-"Burden.pval"
    for(i in 1:n_studies){   #(4)#

        gdata <- subset(study_list[[i]],Group==group & cMAC>min_study_cmac)
        if(gg==1){
        result[[i]]<-gdata
      }else{
        result[[i]]<-rbind(result[[i]],gdata)
      }
        if(nrow(gdata)>0){
        gpvalue0<-unique(unlist(gdata[,pvalues]))
        gpvalue0<-na.omit(gpvalue0)
      }
      gpvalue<-c(gpvalue,gpvalue0)
    }
        burden.cauchy.pval<-CCT(gpvalue)
        cauchy.result1<-data.frame(Group=group,Cauchy.pval=burden.cauchy.pval)
        cauchy.result<-rbind(cauchy.result,cauchy.result1)
        }
    result[["Cauchy.result"]]<- cauchy.result

    return(result)
  }






    # Find transcript group
    transcripts <- lapply(study_list, function(x) unique(x[which(x$Group == group),"Transcript"]) )
    transcripts <- unique(do.call(c, transcripts))
    transcripts <- transcripts[!is.na(transcripts)]

    if(length(transcripts)>0){ #(2)#
      #for(tt in 1:length(transcripts)){ #(3)#
         #print(tt)
         transID <- transcripts[tt]
         n_studies_effective <- NULL
         variant.id <- n.alts <- NULL
         n.site.list<-n.alt.list<-V.list <- U.list <- NULL

      # --------------------------------------------------------------------
      # creating sv.list and V.list
      for(i in 1:n_studies){   #(4)#

          trans_results <- subset(study_list[[i]],Transcript==transID & cMAC>min_study_cmac)
          if(nrow(trans_results)>0){
          n_studies_effective<-c(n_studies_effective,1)
          U.list<-c(U.list,trans_results$Burden.Score)
          V.list<-c(V.list,trans_results$Burden.Variance)
          n.site.list<-c(n.site.list,trans_results$n.site)
          n.alt.list<-c(n.alt.list,trans_results$n.alt)
        }else{
          n_studies_effective<-c(n_studies_effective,0)
          U.list<-c(U.list,NA)
          V.list<-c(V.list,NA)
          n.site.list<-c(n.site.list,NA)
          n.alt.list<-c(n.alt.list,NA)
          }
        }
        U.sum<-sum(U.list,na.rm=T)
        V.sum<-sum(V.list,na.rm=T)

        out <- data.frame(group, transcript=transID, n_studies_effective=sum(n_studies_effective),cMAC=sum(n.alt.list), stringsAsFactors=F)
        colnames(out) <- c("Group", "Transcript","n.studies.contributing", "cMAC")
        nsites<-data.frame(t(n.site.list))
        colnames(nsites)<-paste0("n.sites.",1:length(n.site.list))
        out<-cbind(out,nsites)
        nalt<-data.frame(t(n.alt.list))
        colnames(nalt)<-paste0("n.alts.",1:length(n.alt.list))
        out<-cbind(out,nalt)
        rownames(out) <- group

          if("Burden" %in% test){
            out0<-burden.fun(U.sum=U.sum,V.sum=V.sum)
           out[,names(out0)]<-unlist(out0)
         }
         out$Est.SE<-1/sqrt(out$Burden.Variance)
         out$Est <- out$Burden.Score/out$Burden.Variance

        res <- dplyr::bind_rows(res, out)
      }
    }
  }
  res<-comb_meta_burden_pvalue(data=res,min.cmac=min_meta_cmac)
  return(res)
} #(0)#


####
#### combined p-value from meta-analysis

comb_meta_burden_pvalue<-function(data,pval.col=c("Burden.pval"),min.cmac=10){

result<-list()
pvalues<-names(data)[names(data) %in% pval.col]

data<-subset(data,cMAC>=min.cmac)
#data$Cauchy.pval<-apply(data[,pvalues],1,function(x){x<-na.omit(x);CCT(x)})

genenames<-unique(data$Group)

sum1<-NULL
for (gg in 1:length(genenames)){
genename<-genenames[gg]
gdata<-subset(data,Group==genename &cMAC>=min.cmac)

if(nrow(gdata)>0){
gpvalue<-unique(unlist(gdata[,pvalues]))
gpvalue<-na.omit(gpvalue)
total.cauchy.pval<-CCT(gpvalue)
gdata$Cauchy.anytranscript.pval<-total.cauchy.pval
sum1<-rbind(sum1,gdata)
}
}

return(sum1)
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
