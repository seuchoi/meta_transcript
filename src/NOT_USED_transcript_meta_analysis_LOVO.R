###
### meta-analysis function
### geneid="ENSG00000172967"
rare_variant_meta_analysis_LOVO<-function(study_path_vector=study_path_vector,geneid=NULL, grouping_path_vector=grouping_path_vector,test=c("Burden", "SKAT","SKATO","SMMAT"), min_study_cmac=1,min_meta_cmac=20,rho=c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1),use.anytranscript=T,combine.pval=T){

### output filename
### grouping filename: transcript_id column is required
### any transcript is optional
### available tests

#example=ENSG00000155657;713
#rare_variant_meta_analysis<-function(study_path_vector, test=c("Burden", "SKAT","SKATO","SMMAT"), min_study_cmac=1,rho=c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1),use.anytranscript=T){
#test=c("Burden", "SKAT","SKATO","SMMAT"); min_study_cmac=1;rho=c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1);use.anytranscript=T;combine.pval=T

### geneid is necessary
if(is.null(geneid)){stop("geneid is necessary for leave one variant out analysis")}
###

 n_studies <- length(study_path_vector)
 study_list <- list()
 for(i in c(1:n_studies)){
   load(study_path_vector[i])
   results<-assoc$results[geneid,]
   variantInfo<-assoc$variantInfo[[geneid]]
   covariance_matrix<-assoc$covariance_matrix[[geneid]]
   sub.assoc<-list()
   sub.assoc[["results"]]<-assoc$results[geneid,]
   sub.assoc[["variantInfo"]][[geneid]]<-assoc$variantInfo[[geneid]]
   sub.assoc[["covariance_matrix"]][[geneid]]<-assoc$covariance_matrix[[geneid]]

   study_list[[i]] <- sub.assoc
 }

 groupings <- NULL
 for(i in c(1:n_studies)){
   groupings <- unique(c(groupings, rownames(study_list[[i]]$results)))
 }


transgroup_list<-list()
 for(i in c(1:n_studies)){
   if(use.anytranscript){
     tempgrp0<-get(load(grouping_path_vector[[i]]))
     grp1<-grp0<-subset(tempgrp0,group_id==geneid)
     names(grp1)[which(names(grp1) %in% c("TranscriptID"))]<-names(grp0)[which(names(grp0) %in% c("TranscriptID"))]<-"transcript_id"
     grp1$transcript_id<-grp1$group_id
     #grp0$CANONICAL<-"-"
     grp0 <- grp0 %>% mutate(CANONICAL = "-")
     transgroup_list[[i]]<-rbind(grp0,grp1)
   }else{
    tempgrp0<-get(load(grouping_path_vector[[i]]))
    transgroup_list[[i]]<-subset(tempgrp0,group_id==geneid)
   }
 }


 # Cycle over groupings
 res <- NULL
 num <- 1
 n__groupings <- length(groupings)

#group<-groupings[which(groupings %in% "ENSG00000155657")]
#group<-groupings[1:15]

   # start gene level
   #for(group in groupings){
   group<-groupings
   outnums<-c(1,seq(0,n__groupings,by=10),n__groupings)
   print(num)
   if(num %in% outnums){
     cat('Busy with group', group, 'which is', num, 'out of', n__groupings, '...\n')
   }
   num <- num + 1

      transcripts <- NULL
      for(i in c(1:n_studies)){
        genegroup<-subset(transgroup_list[[i]],group_id==group)
        transcripts <- unique(c(transcripts, genegroup$transcript_id))
        transcripts<-transcripts[!is.na(transcripts)]
        }

if(length(transcripts)>0){
# use variant list
inputvarlist<-c()
for(i in c(1:n_studies)){
varinfodata<-study_list[[i]][["variantInfo"]][[geneid]]
varinfodata<-subset(varinfodata,MAC>0)
varinfodata$variant.id<-paste(varinfodata$chr,varinfodata$pos,varinfodata$ref,varinfodata$alt,sep=":")
inputvarlist<-c(inputvarlist,varinfodata$variant.id)
}
inputvarlist<-unique(inputvarlist)

### annotation variant list
groupingvarlist<-c()
for(i in c(1:n_studies)){
       transgroup_data<-transgroup_list[[i]]
       transgroup_data$variant.id<-paste(transgroup_data$chr,transgroup_data$pos,transgroup_data$ref,transgroup_data$alt,sep=":")
       groupingvarlist<-c(groupingvarlist,transgroup_data$variant.id)
     }
groupingvarlist<-unique(groupingvarlist)

## variants use for LOVO
inputvarlist<-inputvarlist[inputvarlist %in% groupingvarlist]

####
#
for(rm in 1:length(inputvarlist)){

  rmvariant<-inputvarlist[rm]
  translovogroup_list<-list()

  for(i in c(1:n_studies)){
       transgroup_data<-transgroup_list[[i]]
       transgroup_data$variant.id<-paste(transgroup_data$chr,transgroup_data$pos,transgroup_data$ref,transgroup_data$alt,sep=":")
       print(nrow(transgroup_data))
       transgrouprm_data<-subset(transgroup_data,variant.id!=rmvariant)
       print(nrow(transgrouprm_data))
       translovogroup_list[[i]]<-transgrouprm_data
     }

for (tt in 1:length(transcripts)){
  print(tt)
   transID<-transcripts[tt]
   variant.id<-n.alts <- NULL
   n_studies_effective <- NA
   V.list <-sv.list <- list()

   ### needs to fix when one sv list is 0
   for(i in c(1:n_studies)){
          transgroup_data<-translovogroup_list[[i]]
          trans_group_variants<-subset(transgroup_data,transcript_id==transID)
          if(nrow(trans_group_variants)>0){

          trans_group_variants$variant.id<-paste(trans_group_variants$chr,trans_group_variants$pos,trans_group_variants$ref,trans_group_variants$alt,sep=":")

          total_variant<-study_list[[i]]$variantInfo[[group]]
          total_cov<-study_list[[i]]$covariance_matrix[[group]]

          sub_score<-total_variant[total_variant$variant.id %in% trans_group_variants$variant.id,]
          sub_score<-unique(sub_score)
          if(sum(sub_score$MAC)>=min_study_cmac){

          variant.id<-c(variant.id,sub_score$variant.id)
          n.alts<-c(n.alts,sum(sub_score$MAC))

          sub_score<-sub_score[order(sub_score$variant.id),]
          vindex<-sub_score$variant.id
          sv.list[[i]]<-sub_score
          V.list[[i]]<-total_cov[vindex,vindex]
          }
        }
        out<-data.frame(group,transcript=transID,n_studies_effective=length(sv.list),rmvariant=rmvariant,n.site=length(unique(variant.id)),n.alt=sum(n.alts),stringsAsFactors=F)
        colnames(out) <- c("Group", "Transcript","n.studies.contributing","rm.variant", "n.site", "n.alt")
        rownames(out) <- group
        class(out$n.studies.contributing) <- class(out$n.site) <- class(out$n.alt) <- "integer"
        }

        if(length(sv.list)>0){

        effective_studies<-which(!unlist(lapply(sv.list,is.null)))
        n_studies_effective<-length(effective_studies)
      }else{
        effective_studies<-NA
        n_studies_effective<-0
      }
   if(n_studies_effective>0){
     for(i in c(effective_studies)){
       check <- nrow(sv.list[[i]])
         check1 <- sv.list[[i]][,'variant.id'] == colnames(V.list[[i]])
         check2 <- colnames(V.list[[i]]) == rownames(V.list[[i]])
         if(F %in% c(check1, check2)){
             stop("Warning: for cohort ", i, "the variants in the single var file, colnames of Covariance file, or rownames of covariance file do not match.\n")
         }
     }

     variant.list <- NULL
     for(i in c(effective_studies)){
         variant.list <- unique(c(variant.list, sv.list[[i]]$variant.id))
     }
     variant.list <- variant.list[order(variant.list)]
     n.site <- n.variants <- length(variant.list)
     #out$n.site <- n.site
     print(length(variant.list))
     # Reconstruct variant and covariance matrix
     U <- matrix(0, n.variants, 1, dimnames=list(variant.list, 'Score'))
     V <- matrix(0, n.variants, n.variants, dimnames=list(variant.list, variant.list))

     for(i in c(effective_studies)){
       variant.list.cohort <- sv.list[[i]][,'variant.id']
         variant.list.cohort <- variant.list.cohort[order(variant.list.cohort)]
       sv.list[[i]] <- sv.list[[i]][order(sv.list[[i]]$variant.id),]
         U[variant.list.cohort,] <- U[variant.list.cohort,] + sv.list[[i]][,'Score']
       if(is.null(ncol(V.list[[i]]))){
         V[variant.list.cohort, variant.list.cohort] <- matrix(V[variant.list.cohort, variant.list.cohort] + V.list[[i]])
       }else{
           V[variant.list.cohort, variant.list.cohort] <- matrix(V[variant.list.cohort, variant.list.cohort] + V.list[[i]][variant.list.cohort,variant.list.cohort])
       }
     }

     ### Run test
     # Start with Burden, this is also used for SMMAT test, and also for SKAT when implemented later
     U.sum = sum(U[,'Score'])
     V.sum = sum(V)
     GG1 <- rowSums(V)

     # Then using adapted GENESIS script
     ### Run burden test
     burden.pval <- pchisq(U.sum^2/V.sum, df=1, lower.tail=FALSE)
     out[,c("Burden.Score", "Burden.Variance", "Burden.pval")] <- c(U.sum, V.sum, burden.pval)
     class(out$Burden.Score) <- class(out$Burden.Variance) <- class(out$Burden.pval) <- "numeric"

     if("SKAT" %in% test){
       ### Run SKAT
       Q <- sum(U^2)
       SKAT.pval <- NA
       SKAT.pval.method <- NA
       if (mean(abs(V)) >= sqrt(.Machine$double.eps)) {
         pv <- GENESIS:::.regular(Q, V, n.variants)
         SKAT.pval <- pv$pval
         SKAT.pval.method <- pv$method
       }
       out[,c('SKAT.pval', 'SKAT.pval.method')] <- c(SKAT.pval, SKAT.pval.method)
       class(out$SKAT.pval) <- "numeric"
       class(out$SKAT.pval.method) <- "character"
     }

     if("SKATO" %in% test){
       ### Run SKAT
       Q <- sum(U^2)
       SKATO.pval <- NA
       SKATO.pval.method <- NA
       if (mean(abs(V)) >= sqrt(.Machine$double.eps)) {
         res_skato <- GMMAT:::.skato_pval(U = U, V = V, rho = rho, method = "davies")
 		     Burden.Score <- res_skato$Burden.score
 		     Burden.Variance <- res_skato$Burden.var
 		     Burden.pval <- res_skato$Burden.pval
 		     SKAT.pval <- res_skato$SKAT.pval
 		     SKATO.pval <-res_skato$p
 		     SKATO.minp <- res_skato$minp
 		     SKATO.minp.rho <- res_skato$minp.rho
       out[,c("Burden.Score","Burden.Variance","Burden.pval","SKAT.pval","SKATO.pval","SKATO.minp","SKATO.minp.rho")] <- c(Burden.Score,Burden.Variance,Burden.pval,SKAT.pval,SKATO.pval,SKATO.minp,SKATO.minp.rho)
     }}

     if("SMMAT" %in% test){
       ### Run SMMAT
       # Compute burden-adjusted SKAT statistic
       U <- U - GG1*U.sum/V.sum
       Q <- sum(U^2)
       V <- V - tcrossprod(GG1)/V.sum
       # SKAT
       theta.pval <- NA
       theta.pval.method <- NA
                               err <- NA
       if (mean(abs(V)) >= sqrt(.Machine$double.eps)) {
         pv <- GENESIS:::.regular(Q, V, n.variants)
         theta.pval <- pv$pval
         theta.pval.method <- pv$method
         err <- pv$err
       }
       # Fisher's method to combine p-values
       SMMAT.pval <- tryCatch(pchisq(-2*log(burden.pval)-2*log(theta.pval), df=4, lower.tail = FALSE), error = function(e) { NA })
       if(is.na(SMMAT.pval)) {
           err <- 1
         SMMAT.pval <- NA
           SMMAT.pval <- burden.pval
       }
       out[,c("theta.pval", "theta.pval.method", "err", "SMMAT.pval")] <- c(theta.pval, theta.pval.method,  err, SMMAT.pval)
       class(out$theta.pval) <- class(out$err) <- class(out$SMMAT.pval) <- "numeric"
       class(out$theta.pval.method) <- "character"
     }
   }
   res <- dplyr::bind_rows(res, out)
 }
}
}
#return(res)
#}
###combine p-values using cauchy distribution
if(combine.pval){
res<-comb_meta_pvlaue_LOVO(data=res,min.cmac=min_meta_cmac)
}

### add original result
result_gene1<-transcript_meta_analysis_per_gene(study_path_vector=study_path_vector, geneid=geneid,grouping_path_vector=grouping_path_vector,test=c("Burden", "SKAT","SKATO","SMMAT"), min_study_cmac=1,min_meta_cmac=20,rho=c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1),use.anytranscript=T,combine.pval=T)

result_gene1$combined.anytranscript$rm.variant<-"none"
genedata<-result_gene1$combined.anytranscript
genedata<-genedata[,names(res$combined.anytranscript)]

genedatabytrans<-result_gene1$by.transcript
genedatabytrans$rm.variant<-"none"
genedatabytrans<-genedatabytrans[,names(res$by.transcript)]


res$combined.anytranscript<-rbind(genedata,res$combined.anytranscript)
res$by.transcript<-rbind(genedatabytrans,res$by.transcript)


return(res)
}




#data<-data2<-result

comb_meta_pvlaue_LOVO<-function(data,pval.col=c("Burden.pval","SKAT.pval","SKATO.pval","SMMAT.pval"),min.cmac=10){

varlist<-unique(data$rm.variant)
combsubdata<-combanysubdata<-NULL
result<-list()
for (vnum in 1:length(varlist)){

rmvarid<-varlist[vnum]

subdata<-subset(data,rm.variant==rmvarid)
pvalues<-names(subdata)[names(subdata) %in% pval.col]

subdata<-subset(subdata,n.alt>=min.cmac)
subdata$Cauchy.pval<-apply(subdata[,pvalues],1,function(x){x<-na.omit(x);CCT(x)})

genenames<-unique(subdata$Group)

sum1<-NULL
for (gg in 1:length(genenames)){
genename<-genenames[gg]
gdata<-subset(subdata,Group==genename &n.alt>=min.cmac)

if(nrow(gdata)>0){
gpvalue<-unlist(gdata[,pvalues])
gpvalue<-na.omit(gpvalue)
total.cauchy.pval<-CCT(gpvalue)
sum0<-data.frame(Group=genename,Cauchy.anytranscript.pval=total.cauchy.pval)
sum1<-rbind(sum1,sum0)
}
}

gres0<-subset(subdata,Transcript %in% sum1$Group)
gres1<-merge(gres0,sum1,by="Group")

combsubdata<-rbind(combsubdata,subdata)
combanysubdata<-rbind(combanysubdata,gres1)
}

result[["by.transcript"]]<-combsubdata
result[["combined.anytranscript"]]<-combanysubdata

return(result)
}


comb_meta_pvlaue_LOVO_figure<-function(result){


resdata<-result[["combined.anytranscript"]]
rmvar<-unlist(resdata[,"rm.variant"])

staaro<-unlist(resdata[,ncol(resdata)])
pvals<-as.numeric(staaro)


trans_pvals<- -log10(pvals)
range(trans_pvals)
plot(x=c(1:length(trans_pvals)),y=trans_pvals,ylab="-log10(p-value)",xlab="",xlim=c(1,length(trans_pvals)),ylim=c(0,max(trans_pvals+2)),xaxt="n",pch=21,col="black",bg="grey50",cex=1.2)
points(x=1,y=trans_pvals[1],pch=23,col="black",bg="indianred",cex=1.5)
axis(side=1,at=c(1:length(trans_pvals)),labels=rmvar,las=2,cex=0.8)
legend("topleft",legend=c("Original","Leave one variant out"),pch=c(23,21),pt.bg=c("indianred","grey50"),col=c("black","black"),bty="n",cex=1.2)

}
