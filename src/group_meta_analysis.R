library(GENESIS)
library(CompQuadForm)
#library(GMMAT)
library(survey)

vgroupmeta<-function(vgroups,infiles,pvalues=c("Burden.pval","SKAT.pval","SMMAT.pval"),outfile){

for (ff in 1:length(infiles)){

av0<-file.exists(infiles[ff])
print(infiles[ff])
print(av0)

###
if(!any(av0)){stop("Where!")}
###

if(av0){
result<-get(load(infiles[ff]))
result1<-result$by.transcript[,c("Group","Transcript","n.studies.contributing","n.site","n.alt",pvalues)]
names(result1)[3:ncol(result1)]<-paste0(vgroups[ff],".",names(result1)[3:ncol(result1)])
if(ff==1){
print('result 0 defined')
result0<-result1
}else{
result0<-merge(result0,result1,by=c("Group","Transcript"),all=T)
}
}else{result0<-NULL}
}

if(!is.null(result0)){

sum1<-NULL
pvalues<-grep("pval",names(result0))
genenames<-unique(result0$Group)

for (gg in 1:length(genenames)){
genename=genenames[gg]
gdata<-subset(result0,Group==genename)
if(nrow(gdata)>0){
gpvalue<-unlist(gdata[,pvalues])
gpvalue<-na.omit(gpvalue)

#bandaid to larger problem, to discuss w Seung Hoan and Kathy
gpvalue[gpvalue == 1] <- 0.99999999

total.cauchy.pval<-CCT(unique(gpvalue))
sum0<-data.frame(Group=genename,Cauchy.anytranscript.combined.pval=total.cauchy.pval)
sum1<-rbind(sum1,sum0)
}
}

result<-merge(result0,sum1,by="Group")
#change to all groups
save(result,file=outfile)
return(result)
}else{
return(result=NULL)
}
}
