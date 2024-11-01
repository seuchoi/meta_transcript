#####
##### UK biobank

git clone -b dev_choi https://github.com/seuchoi/meta_transcript.git
git pull https://github.com/seuchoi/meta_transcript.git


cp /mnt/project/schoi/rpackages/rpackages4_2_2.tar.gz ./

tar -xvf rpackages4_2_2.tar.gz


.libPaths(c("rpackages4_2_2",.libPaths()))
R.utils::sourceDirectory("/opt/notebooks/meta_transcript/src/")
install.packages(c("CompQuadForm","survey"))
library(GENESIS)
library(CompQuadForm)
library(survey)

###
### hclof_noflag
for (num in 1:22){

study_path<-paste0("/mnt/project/schoi/association/hclof_noflag_missense0.8/CAD_HARD_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr",num,".RData")
grouping_path<-paste0("/mnt/project/schoi/annotation/group/ukb23156_c",num,"_genotype_variant_QCed_merged.annotated.vep105.gz.hclof_noflag_gnomad_POPMAX0.001.RData")
outfile<-paste0("CAD_HARD_Burden_hclof_noflag_POPMAX0.001_chr",num,".RData")
result<-transcript_single_analysis_nocovout(study_path=study_path,
                                     grouping_path=grouping_path,test=c("Burden"),
                                     min_study_cmac=1,min_meta_cmac=20,
                                     use.anytranscript=T,combine.pval=T)

save(result,file=outfile)
}

dx mkdir schoi/association/burden/hclof_noflag
dx upload CAD_HARD_Burden_hclof_noflag_POPMAX0.001*.RData --path cad_rvas:/schoi/association/burden/hclof_noflag/
###
### hclof_noflag_missense0.8

for (num in 1:22){

study_path<-paste0("/mnt/project/schoi/association/hclof_noflag_missense0.8/CAD_HARD_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr",num,".RData")
grouping_path<-paste0("/mnt/project/schoi/annotation/group/ukb23156_c",num,"_genotype_variant_QCed_merged.annotated.vep105.gz.hclof_noflag_missense0.8_gnomad_POPMAX0.001.RData")
outfile<-paste0("CAD_HARD_Burden_hclof_noflag_missense0.8_POPMAX0.001_chr",num,".RData")
result<-transcript_single_analysis_nocovout(study_path=study_path,
                                     grouping_path=grouping_path,test=c("Burden"),
                                     min_study_cmac=1,min_meta_cmac=20,
                                     use.anytranscript=T,combine.pval=T)

save(result,file=outfile)
}
dx mkdir schoi/association/burden/hclof_noflag_missense0.8
dx upload CAD_HARD_Burden_hclof_noflag_missense0.8_POPMAX0.001*.RData --path cad_rvas:/schoi/association/burden/hclof_noflag_missense0.8/


########
######## migen
R.utils::sourceDirectory("/medpop/esp2/schoi/software/meta_transcript/src")
library(GENESIS)
library(CompQuadForm)
library(survey)


setwd("/medpop/esp2/schoi/migen/migen/result/association/burden")
for (num in 1:22){

study_path<-paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/MIGEN_V13.hclof_noflag_missense0.8_POPMAX0.001.chr",num,".RData")
grouping_path<-paste0("/medpop/esp2/skoyama/passing/migen/ws01/70_annotation/out/grouping/hclof_noflag_POPMAX0.001/MIGEN_V13_sampleQCed_chr",num,".pvar.vep105.gz.hclof_noflag_POPMAX0.001.RData")
outfile<-paste0("MIGEN_V13_Burden_hclof_noflag_POPMAX0.001_chr",num,".RData")
result<-transcript_single_analysis_nocovout(study_path=study_path,
                                     grouping_path=grouping_path,test=c("Burden"),
                                     min_study_cmac=1,min_meta_cmac=20,
                                     use.anytranscript=T,combine.pval=T)

save(result,file=outfile)
}


setwd("/medpop/esp2/schoi/migen/migen/result/association/burden")
for (num in 1:22){

study_path<-paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/MIGEN_V13.hclof_noflag_missense0.8_POPMAX0.001.chr",num,".RData")
grouping_path<-paste0("/medpop/esp2/skoyama/passing/migen/ws01/70_annotation/out/grouping/hclof_noflag_missense0.8_POPMAX0.001/MIGEN_V13_sampleQCed_chr",num,".pvar.vep105.gz.hclof_noflag_missense0.8_POPMAX0.001.RData")
outfile<-paste0("MIGEN_V13_Burden_hclof_noflag_missense0.8_POPMAX0.001_chr",num,".RData")
result<-transcript_single_analysis_nocovout(study_path=study_path,
                                     grouping_path=grouping_path,test=c("Burden"),
                                     min_study_cmac=1,min_meta_cmac=20,
                                     use.anytranscript=T,combine.pval=T)

save(result,file=outfile)
}


#####
##### AoU
.libPaths(c("rpackages4_1_3",.libPaths()))
R.utils::sourceDirectory("/home/jupyter/workspaces/cadrvas/meta_transcript/src")
library(GENESIS)
library(CompQuadForm)
library(survey)


setwd("/home/jupyter/workspaces/cadrvas/result/burden/hclof_noflag")
for (num in 1:22){

study_path<-paste0("/home/jupyter/workspaces/cadrvas/result/cad_185k_CAD_chr",num,".RData")
grouping_path<-paste0("/home/jupyter/workspaces/cadrvas/annotation/AoU_250K_exome_annotation_VEP105_chr",num,".vcf.gz.hclof_noflag_POPMAX0.001.RData")
outfile<-paste0("AoU_Burden_hclof_noflag_POPMAX0.001_chr",num,".RData")
result<-transcript_single_analysis_nocovout(study_path=study_path,
                                     grouping_path=grouping_path,test=c("Burden"),
                                     min_study_cmac=1,min_meta_cmac=20,
                                     use.anytranscript=T,combine.pval=T)
message(paste0("chr",num," is done"))
save(result,file=outfile)
}




setwd("/home/jupyter/workspaces/cadrvas/result/burden/hclof_noflag_missense0.8")
for (num in 1:22){

study_path<-paste0("/home/jupyter/workspaces/cadrvas/result/cad_185k_CAD_chr",num,".RData")
grouping_path<-paste0("/home/jupyter/workspaces/cadrvas/annotation/AoU_250K_exome_annotation_VEP105_chr",num,".vcf.gz.hclof_noflag_missense0.8_7tools_POPMAX0.001.RData")
outfile<-paste0("AoU_Burden_hclof_noflag_missense0.8_POPMAX0.001_chr",num,".RData")
result<-transcript_single_analysis_nocovout(study_path=study_path,
                                     grouping_path=grouping_path,test=c("Burden"),
                                     min_study_cmac=1,min_meta_cmac=20,
                                     use.anytranscript=T,combine.pval=T)
message(paste0("chr",num," is done"))
save(result,file=outfile)
}


tar -czvf test.tar.gz /home/jupyter/workspaces/cadrvas/result/burden


####
#### UK biobank
cd /medpop/afib/

use .python-3.8.3
source dxpy/bin/activate

#dx login
dx select "cad_rvas"

dx download cad_rvas:/schoi/association/burden/hclof_noflag/* ./
dx download cad_rvas:/schoi/association/burden/hclof_noflag_missense0.8/*


dx download cad_rvas:/schoi/association/burden/hclof_noflag_missense0.8/UKBB_migen_EUR_MI_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr2_SPA.RData 


######
###### cohort meta-analysis

R.utils::sourceDirectory("/medpop/esp2/schoi/software/meta_transcript/src")
library(GENESIS)
library(CompQuadForm)
library(survey)

setwd("/medpop/esp2/schoi/migen/meta/hclof_noflag")
for (num in 1:22){
ukbb<-paste0("/medpop/esp2/schoi/migen/ukbb/result/burden/hclof_noflag/CAD_HARD_Burden_hclof_noflag_POPMAX0.001_chr",num,".RData")
aou<-paste0("/medpop/esp2/schoi/migen/aou/hclof_noflag/AoU_Burden_hclof_noflag_POPMAX0.001_chr",num,".RData")
migen<-paste0("/medpop/esp2/schoi/migen/migen/result/association/burden/hclof_noflag/MIGEN_V13_Burden_hclof_noflag_POPMAX0.001_chr",num,".RData")
outfile<-paste0("UKBB_AoU_MIGEN_Burden_hclof_noflag_POPMAX0.001_chr",num,".RData")
study_path<-c(ukbb,aou,migen)

result<-transcript_single_analysis_meta_nocovout(study_path=study_path,test=c("Burden"),
                                   min_study_cmac=20,min_meta_cmac=20,
                                   use.anytranscript=T)

save(result,file=outfile)
}


setwd("/medpop/esp2/schoi/migen/meta/hclof_noflag_missense0.8")
for (num in 1:22){
ukbb<-paste0("/medpop/esp2/schoi/migen/ukbb/result/burden/hclof_noflag_missense0.8/CAD_HARD_Burden_hclof_noflag_missense0.8_POPMAX0.001_chr",num,".RData")
aou<-paste0("/medpop/esp2/schoi/migen/aou/hclof_noflag_missense0.8/AoU_Burden_hclof_noflag_missense0.8_POPMAX0.001_chr",num,".RData")
migen<-paste0("/medpop/esp2/schoi/migen/migen/result/association/burden/hclof_noflag_missense0.8/MIGEN_V13_Burden_hclof_noflag_missense0.8_POPMAX0.001_chr",num,".RData")
outfile<-paste0("UKBB_AoU_MIGEN_Burden_hclof_noflag_missense0.8_POPMAX0.001_chr",num,".RData")
study_path<-c(ukbb,aou,migen)

result<-transcript_single_analysis_meta_nocovout(study_path=study_path,test=c("Burden"),
                                   min_study_cmac=20,min_meta_cmac=20,
                                   use.anytranscript=T)

save(result,file=outfile)
}

####
#### check the number of genes
all<-NULL
for (ii in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/meta/hclof_noflag/UKBB_AoU_MIGEN_Burden_hclof_noflag_POPMAX0.001_chr",ii,".RData")
aa<-load(filename)

all<-rbind(all,result)
}
length(unique(all$Group))






######
######


for (num in 1:22){


study_path<-c(paste0("/medpop/esp2/schoi/migen/meta/hclof_noflag/UKBB_AoU_MIGEN_Burden_hclof_noflag_POPMAX0.001_chr",num,".RData"),paste0("/medpop/esp2/schoi/migen/meta/hclof_noflag_missense0.8/UKBB_AoU_MIGEN_Burden_hclof_noflag_missense0.8_POPMAX0.001_chr",num,".RData"))

result<-transcript_cross_annot_meta(study_path=study_path,test=c("Burden"),
                                   min_study_cmac=1,min_meta_cmac=20,
                                   use.anytranscript=T)

outfile<-paste0("/medpop/esp2/schoi/migen/meta/UKBB_AoU_MIGEN_Burden_Cauchy_chr",num,".RData")
save(result,file=outfile)
}

######
###### manhattan plot
allres<-NULL
for (num in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/meta/UKBB_AoU_MIGEN_Burden_Cauchy_chr",num,".RData")
load(filename)
sum0<-result[["Cauchy.result"]]
sum0$chr<-num
allres<-rbind(allres,sum0)
}

####
#### manhattan plot
aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(allres,gencode[,c("geneid","gene_name","start")],by.x="Group",by.y="geneid")
siggenes<-subset(comb1,Cauchy.pval<3.189589e-06)
siggenelist<-siggenes$gene_name
setwd("/medpop/esp2/schoi/migen/meta/figure")

library(qqman)

png("UKBB_AoU_MIGEN_Burden_Cauchy_manhattan.png",width = 1000, height = 500,res=110)
manhattan(comb1,chr="chr",bp="start",p="Cauchy.pval",snp="gene_name",col = c("#1b9e77","#d95f02"),genomewideline = -log10(3.189589e-06),
    suggestiveline = -log10(6.379178e-06),annotatePval=3.189589e-06,annotateTop = FALSE)
dev.off()


png("UKBB_AoU_MIGEN_Burden_Cauchy_qq.png",width = 1000, height = 1000,res=110)
qq(comb1$Cauchy.pval)
alpha<-median(qchisq(1-comb1$Cauchy.pval,1))/qchisq(0.5,1)
text(0.5,12, paste("lambda","=",  signif(alpha, digits = 3)) )
dev.off()



### test1
all<-NULL
for (ii in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/ukbb/result/burden/hclof_noflag_missense0.8/CAD_HARD_Burden_hclof_noflag_missense0.8_POPMAX0.001_chr",num,".RData")
aa<-load(filename)
all<-rbind(all,result$combined.anytranscript)
}
nrow(all)



}



####
#### SPA
cd /medpop/afib/

use .python-3.8.3
source dxpy/bin/activate

#dx login
dx select "cad_rvas"

### download again!
### CAD
cd /medpop/esp2/schoi/migen/ukbb/result/burden/spa/cad/hclof_noflag_POPMAX0.001
dx download -f cad_rvas:/schoi/association/spa/hclof_noflag/UKBB_migen_*_CAD_hclof_noflag_POPMAX0.001_chr*_SPA.RData
dx ls cad_rvas:/schoi/association/spa/hclof_noflag_missense0.8/UKBB_migen_*_CAD_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr*_SPA.RData
cd /medpop/esp2/schoi/migen/ukbb/result/burden/spa/cad/hclof_noflag_missense0.8_7tools_POPMAX0.001
dx download -f cad_rvas:/schoi/association/spa/hclof_noflag_missense0.8/UKBB_migen_*_CAD_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr*_SPA.RData

### MI
cd /medpop/esp2/schoi/migen/ukbb/result/burden/spa/mi/hclof_noflag_POPMAX0.001
dx ls cad_rvas:/schoi/association/spa/hclof_noflag/UKBB_migen_*_MI_hclof_noflag_POPMAX0.001_chr*_SPA.RData
dx download -f cad_rvas:/schoi/association/spa/hclof_noflag/UKBB_migen_*_MI_hclof_noflag_POPMAX0.001_chr*_SPA.RData

cd /medpop/esp2/schoi/migen/ukbb/result/burden/spa/mi/hclof_noflag_missense0.8_7tools_POPMAX0.001
dx ls cad_rvas:/schoi/association/spa/hclof_noflag_missense0.8/UKBB_migen_*_MI_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr*_SPA.RData
dx download -f cad_rvas:/schoi/association/spa/hclof_noflag_missense0.8/UKBB_migen_*_MI_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr*_SPA.RData

#### EOMI
cd /medpop/esp2/schoi/migen/ukbb/result/burden/spa/eomi/hclof_noflag_POPMAX0.001
dx ls cad_rvas:/schoi/association/spa/hclof_noflag/UKBB_migen_*_EOMI_hclof_noflag_POPMAX0.001_chr*_SPA.RData
dx download -f cad_rvas:/schoi/association/spa/hclof_noflag/UKBB_migen_*_EOMI_hclof_noflag_POPMAX0.001_chr*_SPA.RData

cd /medpop/esp2/schoi/migen/ukbb/result/burden/spa/eomi/hclof_noflag_missense0.8_7tools_POPMAX0.001
dx ls cad_rvas:/schoi/association/spa/hclof_noflag_missense0.8/UKBB_migen_*_EOMI_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr*_SPA.RData
dx download -f cad_rvas:/schoi/association/spa/hclof_noflag_missense0.8/UKBB_migen_*_EOMI_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr*_SPA.RData



#### All of us
#### SPA
R.utils::sourceDirectory("/medpop/esp2/schoi/software/meta_transcript/src")

#### prepare
setwd("/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_POPMAX0.001")

traits<-c("cad","mi","eomi")
for (tt in 1:length(traits)){
trait<-traits[tt]    

load(paste0(trait,"_hclof_noflag_POPMAX0.001_chrAll_SPA.RData"))
res0<-assoc
for (num in 1:22){
assoc<-list()
results<-subset(res0,chr==num)
assoc[["results"]]<-results
save(assoc,file=paste0(trait,"_hclof_noflag_POPMAX0.001_chr",num,"_SPA.RData"))
}
}


#### prepare
setwd("/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_missense0.8_POPMAX0.001")

traits<-c("cad","mi","eomi")
for (tt in 1:length(traits)){
trait<-traits[tt]    

load(paste0(trait,"_hclof_noflag_missense0.8_POPMAX0.001_chrAll_SPA.RData"))
res0<-assoc
for (num in 1:22){
assoc<-list()
results<-subset(res0,chr==num)
assoc[["results"]]<-results
save(assoc,file=paste0(trait,"_hclof_noflag_missense0.8_POPMAX0.001_chr",num,"_SPA.RData"))
}
}





args=(commandArgs(TRUE))
input1=as.character(args[1])
input2=as.character(args[2])
input3=as.character(args[3])
input4=as.character(args[4])
output=as.character(args[5])
source("/medpop/esp2/schoi/software/meta_transcript/src/rarevar_tests.R")
source("/medpop/esp2/schoi/software/meta_transcript/src/transcript_single_analysis_meta_nocov_spa.R")
library(GENESIS)
library(CompQuadForm)
library(survey)

study_path<-c(input1,input2,input3,input4)

result<-transcript_single_analysis_meta_nocovout_spa(study_path=study_path,test=c("Burden"),min_study_cmac=20,min_meta_cmac=20)
save(result,file=output)

sessionInfo()
quit("no")


#!/bin/bash -l
#$ -l h_rt=4:00:00
#$ -l h_vmem=16g
#$ -N transcript_meta_analysis
#$ -cwd
#$ -j y
#$ -o log/$JOB_NAME.$JOB_ID.log
echo "=========================================================="
echo "Starting on       : $(date)"
echo "Running on node   : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID    : $JOB_ID"
echo "Current job name  : $JOB_NAME"
echo "Task index number : chrom ${num}"
echo "=========================================================="

reuse .perl-5.28.0
reuse .tabix-0.2.6
reuse .samtools-1.8
reuse .openssl-1.0.2g
reuse .jq-1.5
reuse R-4.0

###
cd /medpop/esp2/schoi/migen/meta/script/
output2="$(basename ${output})"

R CMD BATCH "--args ${input1} ${input2} ${input3} ${input4} ${output}" /medpop/esp2/schoi/migen/meta/script/transcript_meta_analysis.R out/${output2}.out

transcript_single_analysis_meta_nocov_spa.R
echo "=========================================================="
echo "Finished on       : $(date)"
echo "Task index number : chrom ${num}"
echo "=========================================================="


cd /medpop/esp2/schoi/migen/meta/script


#####hclof_noflag_POPMAX0.001
#### CAD
for num in {1..22}
do
in1=/medpop/esp2/schoi/migen/ukbb/result/burden/spa/cad/hclof_noflag_POPMAX0.001/UKBB_migen_ALL_CAD_hclof_noflag_POPMAX0.001_chr${num}_SPA.RData
in2=/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_POPMAX0.001/cad_hclof_noflag_POPMAX0.001_chr${num}_SPA.RData
in3=/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_POPMAX0.001/spa/MIGEN_V13.CHD.hclof_noflag_POPMAX0.001_chr${num}_SPA.RData
in4=/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_POPMAX0.001/MGBB.cad.hclof_noflag_POPMAX0.001.chr${num}_SPA.RData
out1=/medpop/esp2/schoi/migen/meta/hclof_noflag/cad_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_POPMAX0.001_SPA_chr${num}.RData

qsub -N chr${num}_hclof_noflag_meta -v input1=${in1},input2=${in2},input3=${in3},input4=${in4},output=${out1} transcript_meta_analysis.qsub

done

### MI

for num in {1..22}
do
in1=/medpop/esp2/schoi/migen/ukbb/result/burden/spa/mi/hclof_noflag_POPMAX0.001/UKBB_migen_ALL_MI_hclof_noflag_POPMAX0.001_chr${num}_SPA.RData
in2=/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_POPMAX0.001/mi_hclof_noflag_POPMAX0.001_chr${num}_SPA.RData
in3=/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_POPMAX0.001/spa/MIGEN_V13.MI.hclof_noflag_POPMAX0.001_chr${num}_SPA.RData
in4=/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_POPMAX0.001/MGBB.mi.hclof_noflag_POPMAX0.001.chr${num}_SPA.RData
out1=/medpop/esp2/schoi/migen/meta/hclof_noflag/mi_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_POPMAX0.001_SPA_chr${num}.RData

qsub -N chr${num}_mi_hclof_noflag_meta -v input1=${in1},input2=${in2},input3=${in3},input4=${in4},output=${out1} transcript_meta_analysis.qsub

done

### eomi
for num in {1..22}
do
in1=/medpop/esp2/schoi/migen/ukbb/result/burden/spa/eomi/hclof_noflag_POPMAX0.001/UKBB_migen_ALL_EOMI_hclof_noflag_POPMAX0.001_chr${num}_SPA.RData
in2=/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_POPMAX0.001/eomi_hclof_noflag_POPMAX0.001_chr${num}_SPA.RData
in3=/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_POPMAX0.001/spa/MIGEN_V13.EOMI.hclof_noflag_POPMAX0.001_chr${num}_SPA.RData
in4=/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_POPMAX0.001/MGBB.eomi.hclof_noflag_POPMAX0.001.chr${num}_SPA.RData
out1=/medpop/esp2/schoi/migen/meta/hclof_noflag/eomi_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_POPMAX0.001_SPA_chr${num}.RData

qsub -N chr${num}_eomi_hclof_noflag_meta -v input1=${in1},input2=${in2},input3=${in3},input4=${in4},output=${out1} transcript_meta_analysis.qsub

done

#####hclof_noflag_missense0.8_POPMAX0.001
#### CAD
for num in {1..22}
do
in1=/medpop/esp2/schoi/migen/ukbb/result/burden/spa/cad/hclof_noflag_missense0.8_7tools_POPMAX0.001/UKBB_migen_ALL_CAD_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr${num}_SPA.RData
in2=/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_missense0.8_POPMAX0.001/cad_hclof_noflag_missense0.8_POPMAX0.001_chr${num}_SPA.RData
in3=/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/spa/MIGEN_V13.CHD.hclof_noflag_missense0.8_POPMAX0.001.chr${num}_SPA.RData
in4=/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_missense0.8_POPMAX0.001/MGBB.cad.hclof_noflag_missense0.8_POPMAX0.001.chr${num}_SPA.RData
out1=/medpop/esp2/schoi/migen/meta/hclof_noflag_missense0.8/UKBB_AoU_MIGEN_MGBB/cad/cad_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_missense0.8_POPMAX0.001_SPA_chr${num}.RData

qsub -N chr${num}_hclof_noflag_missense0.8_meta -v input1=${in1},input2=${in2},input3=${in3},input4=${in4},output=${out1} transcript_meta_analysis.qsub

done

### MI
for num in {1..22}
do
in1=/medpop/esp2/schoi/migen/ukbb/result/burden/spa/mi/hclof_noflag_missense0.8_7tools_POPMAX0.001/UKBB_migen_ALL_MI_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr${num}_SPA.RData
in2=/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_missense0.8_POPMAX0.001/mi_hclof_noflag_missense0.8_POPMAX0.001_chr${num}_SPA.RData
in3=/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/spa/MIGEN_V13.MI.hclof_noflag_missense0.8_POPMAX0.001.chr${num}_SPA.RData
in4=/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_missense0.8_POPMAX0.001/MGBB.mi.hclof_noflag_missense0.8_POPMAX0.001.chr${num}_SPA.RData
ls ${in1};
ls ${in2};
ls ${in3};
ls ${in4};
done

for num in {1..22}
do
in1=/medpop/esp2/schoi/migen/ukbb/result/burden/spa/mi/hclof_noflag_missense0.8_7tools_POPMAX0.001/UKBB_migen_ALL_MI_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr${num}_SPA.RData
in2=/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_missense0.8_POPMAX0.001/mi_hclof_noflag_missense0.8_POPMAX0.001_chr${num}_SPA.RData
in3=/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/spa/MIGEN_V13.MI.hclof_noflag_missense0.8_POPMAX0.001.chr${num}_SPA.RData
in4=/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_missense0.8_POPMAX0.001/MGBB.mi.hclof_noflag_missense0.8_POPMAX0.001.chr${num}_SPA.RData
out1=/medpop/esp2/schoi/migen/meta/hclof_noflag_missense0.8/UKBB_AoU_MIGEN_MGBB/mi/mi_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_missense0.8_POPMAX0.001_SPA_chr${num}.RData

qsub -N chr${num}_mi_hclof_noflag_missense0.8_meta -v input1=${in1},input2=${in2},input3=${in3},input4=${in4},output=${out1} transcript_meta_analysis.qsub

done

### eomi
for num in {1..22}
do
in1=/medpop/esp2/schoi/migen/ukbb/result/burden/spa/eomi/hclof_noflag_missense0.8_7tools_POPMAX0.001/UKBB_migen_ALL_EOMI_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr${num}_SPA.RData
in2=/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_missense0.8_POPMAX0.001/eomi_hclof_noflag_missense0.8_POPMAX0.001_chr${num}_SPA.RData
in3=/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/spa/MIGEN_V13.EOMI.hclof_noflag_missense0.8_POPMAX0.001.chr${num}_SPA.RData
in4=/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_missense0.8_POPMAX0.001/MGBB.eomi.hclof_noflag_missense0.8_POPMAX0.001.chr${num}_SPA.RData
out1=/medpop/esp2/schoi/migen/meta/hclof_noflag_missense0.8/eomi_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_missense0.8_POPMAX0.001_SPA_chr${num}.RData


qsub -N chr${num}_eomi_hclof_noflag_missense0.8_meta -v input1=${in1},input2=${in2},input3=${in3},input4=${in4},output=${out1} transcript_meta_analysis.qsub

done


## MI

for num in {1..22}
do
in1=/medpop/esp2/schoi/migen/ukbb/result/burden/spa/mi/hclof_noflag_missense0.8_7tools_POPMAX0.001/UKBB_migen_ALL_MI_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr${num}_SPA.RData
in2=/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_missense0.8_POPMAX0.001/mi_hclof_noflag_missense0.8_POPMAX0.001_chr${num}_SPA.RData
in3=/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/spa/MIGEN_V13.MI.hclof_noflag_missense0.8_POPMAX0.001.chr${num}.RData
in4=/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_missense0.8_POPMAX0.001/MGBB.mi.hclof_noflag_missense0.8_POPMAX0.001.chr${num}_SPA.RData
out1=/medpop/esp2/schoi/migen/meta/hclof_noflag_missense0.8/UKBB_AoU_MIGEN_MGBB/mi/mi_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_missense0.8_POPMAX0.001_SPA_chr${num}.RData

qsub -N chr${num}_mi_hclof_noflag_missense0.8_meta -v input1=${in1},input2=${in2},input3=${in3},input4=${in4},output=${out1} transcript_meta_analysis.qsub

done




### mi example
### hclof_noflag_missense0.8
setwd("/medpop/esp2/schoi/migen/meta/hclof_noflag_missense0.8/UKBB_AoU_MIGEN_MGBB/mi")
for (num in 1:22){
ukbb<-paste0("/medpop/esp2/schoi/migen/ukbb/result/burden/spa/mi/hclof_noflag_missense0.8_7tools_POPMAX0.001/UKBB_migen_ALL_MI_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr",num,"_SPA.RData")
aou<-paste0("/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_missense0.8_POPMAX0.001/mi_hclof_noflag_missense0.8_POPMAX0.001_chr",num,"_SPA.RData")
migen<-paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/MIGEN_V13.MI.hclof_noflag_missense0.8_POPMAX0.001.chr",num,".RData")
mgbb<-paste0("/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_missense0.8_POPMAX0.001/MGBB.mi.hclof_noflag_missense0.8_POPMAX0.001.chr",num,"_SPA.RData")
outfile<-paste0("/medpop/esp2/schoi/migen/meta/hclof_noflag_missense0.8/UKBB_AoU_MIGEN_MGBB/mi/mi_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_missense0.8_POPMAX0.001_SPA_chr",num,".RData")

study_path<-c(ukbb,aou,migen,mgbb)

result<-transcript_single_analysis_meta_nocovout_spa(study_path=study_path,test=c("Burden"),
                                   min_study_cmac=20,min_meta_cmac=20,
                                   use.anytranscript=T)

save(result,file=outfile)
}





###CAD
source("/medpop/esp2/schoi/software/meta_transcript/src/rarevar_tests.R")
source("/medpop/esp2/schoi/software/meta_transcript/src/transcript_cross_annot_meta.R")

for (num in 1:22){

study_path<-c(paste0("/medpop/esp2/schoi/migen/meta/hclof_noflag/UKBB_AoU_MIGEN_MGBB/cad/cad_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_POPMAX0.001_SPA_chr",num,".RData"),paste0("/medpop/esp2/schoi/migen/meta/hclof_noflag_missense0.8/UKBB_AoU_MIGEN_MGBB/cad/cad_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_missense0.8_POPMAX0.001_SPA_chr",num,".RData"))

result<-transcript_cross_annot_meta(study_path=study_path,test=c("Burden"),
                                   min_study_cmac=20,min_meta_cmac=20,
                                   use.anytranscript=T)

outfile<-paste0("/medpop/esp2/schoi/migen/meta/cad_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",num,".RData")
save(result,file=outfile)
}

########
########
R.utils::sourceDirectory("/medpop/esp2/schoi/software/meta_transcript/src")

library(GENESIS)
library(CompQuadFrom)
library(survey)

######
###### manhattan plot
allres<-NULL
for (num in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/meta/cad_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",num,".RData")
load(filename)
sum0<-result[["Cauchy.result"]]
sum0$chr<-num
allres<-rbind(allres,sum0)
}
nrow(allres)

####
#### manhattan plot
aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(allres,gencode[,c("geneid","gene_name","start")],by.x="Group",by.y="geneid")
gthreshold<-0.05/nrow(comb1)
sthreshold<-0.1/nrow(comb1)

siggenes<-subset(comb1,Cauchy.pval<gthreshold)
siggenelist<-siggenes$gene_name

setwd("/medpop/esp2/schoi/migen/meta/figure")
library(qqman)

png("UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_manhattan.png",width = 1000, height = 500,res=110)
manhattan(comb1,chr="chr",bp="start",p="Cauchy.pval",snp="gene_name",col = c("#1b9e77","#d95f02"),genomewideline = -log10(gthreshold),
    suggestiveline = -log10(sthreshold),annotatePval=gthreshold,annotateTop = FALSE)
dev.off()


png("UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_qq.png",width = 1000, height = 1000,res=110)
qq(comb1$Cauchy.pval)
alpha<-median(qchisq(1-comb1$Cauchy.pval,1))/qchisq(0.5,1)
text(0.5,12, paste("lambda","=",  signif(alpha, digits = 3)) )
dev.off()


###### LOF only
allres<-NULL
for (num in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/meta/cad_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",num,".RData")
load(filename)
sum0<-result[["mask1"]]
sum0$chr<-num
#sum1<-subset(sum0,Group %in% siggenes$Group)
allres<-rbind(allres,sum0)
}
allres$type<-"hclof_noflag"
lof<-allres


####
#### manhattan plot
aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(lof,gencode[,c("geneid","gene_name","start")],by.x="Group",by.y="geneid")

setwd("/medpop/esp2/schoi/migen/meta/figure")

library(qqman)

comb2<-unique(comb1[,c("Group","chr","start","gene_name","Cauchy.anytranscript.pval")])

png("UKBB_AoU_MIGEN_MGBB_hclof_noflag_Burden_Cauchy_manhattan.png",width = 1000, height = 500,res=110)
manhattan(comb2,chr="chr",bp="start",p="Cauchy.anytranscript.pval",snp="gene_name",col = c("#1b9e77","#d95f02"),genomewideline = -log10(gthreshold),
    suggestiveline = -log10(sthreshold),annotatePval=gthreshold,annotateTop = FALSE)
dev.off()


png("UKBB_AoU_MIGEN_MGBB_hclof_noflag_Burden_Cauchy_qq.png",width = 1000, height = 1000,res=110)
qq(comb2$Cauchy.anytranscript.pval)
alpha<-median(qchisq(1-comb2$Cauchy.anytranscript.pval,1))/qchisq(0.5,1)
text(0.5,12, paste("lambda","=",  signif(alpha, digits = 3)) )
dev.off()

loflist<-subset(comb2,Cauchy.anytranscript.pval<gthreshold)


######
###### hclof_noflag_missense0.8
allres<-NULL
for (num in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/meta/cad_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",num,".RData")
load(filename)
sum0<-result[["mask2"]]
sum0$chr<-num
#sum1<-subset(sum0,Group %in% siggenes$Group)
allres<-rbind(allres,sum0)
}
allres$type<-"hclof_noflag_missense0.8"
lof_missense<-allres

####
####
aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(lof_missense,gencode[,c("geneid","gene_name","start")],by.x="Group",by.y="geneid")

setwd("/medpop/esp2/schoi/migen/meta/figure")

library(qqman)

comb2<-unique(comb1[,c("Group","chr","start","gene_name","Cauchy.anytranscript.pval")])

png("UKBB_AoU_MIGEN_MGBB_hclof_noflag_missense0.8_Burden_Cauchy_manhattan.png",width = 1000, height = 500,res=110)
manhattan(comb2,chr="chr",bp="start",p="Cauchy.anytranscript.pval",snp="gene_name",col = c("#1b9e77","#d95f02"),genomewideline = -log10(gthreshold),
    suggestiveline = -log10(sthreshold),annotatePval=gthreshold,annotateTop = FALSE)
dev.off()


png("UKBB_AoU_MIGEN_MGBB_hclof_noflag_missense0.8_Burden_Cauchy_qq.png",width = 1000, height = 1000,res=110)
qq(comb2$Cauchy.anytranscript.pval)
alpha<-median(qchisq(1-comb2$Cauchy.anytranscript.pval,1))/qchisq(0.5,1)
text(0.5,12, paste("lambda","=",  signif(alpha, digits = 3)) )
dev.off()

####



#######
###eomi
source("/medpop/esp2/schoi/software/meta_transcript/src/rarevar_tests.R")
source("/medpop/esp2/schoi/software/meta_transcript/src/transcript_cross_annot_meta.R")

for (num in 1:22){

study_path<-c(paste0("/medpop/esp2/schoi/migen/meta/hclof_noflag/UKBB_AoU_MIGEN_MGBB/eomi/eomi_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_POPMAX0.001_SPA_chr",num,".RData"),paste0("/medpop/esp2/schoi/migen/meta/hclof_noflag_missense0.8/UKBB_AoU_MIGEN_MGBB/eomi/eomi_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_missense0.8_POPMAX0.001_SPA_chr",num,".RData"))

result<-transcript_cross_annot_meta(study_path=study_path,test=c("Burden"),
                                   min_study_cmac=20,min_meta_cmac=20,
                                   use.anytranscript=T)

outfile<-paste0("/medpop/esp2/schoi/migen/meta/eomi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",num,".RData")
save(result,file=outfile)
}

########
########
R.utils::sourceDirectory("/medpop/esp2/schoi/software/meta_transcript/src")

library(GENESIS)
library(CompQuadFrom)
library(survey)

######
###### manhattan plot
allres<-NULL
for (num in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/meta/eomi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",num,".RData")
load(filename)
sum0<-result[["Cauchy.result"]]
sum0$chr<-num
allres<-rbind(allres,sum0)
}
nrow(allres)

####
#### manhattan plot
aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(allres,gencode[,c("geneid","gene_name","start")],by.x="Group",by.y="geneid")
gthreshold<-0.05/nrow(comb1)
sthreshold<-0.1/nrow(comb1)

siggenes<-subset(comb1,Cauchy.pval<gthreshold)
siggenelist<-siggenes$gene_name

setwd("/medpop/esp2/schoi/migen/meta/figure")
library(qqman)

png("eomi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_manhattan.png",width = 1000, height = 500,res=110)
manhattan(comb1,chr="chr",bp="start",p="Cauchy.pval",snp="gene_name",col = c("#1b9e77","#d95f02"),genomewideline = -log10(gthreshold),
    suggestiveline = -log10(sthreshold),annotatePval=gthreshold,annotateTop = FALSE)
dev.off()


png("eomi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_qq.png",width = 1000, height = 1000,res=110)
qq(comb1$Cauchy.pval)
alpha<-median(qchisq(1-comb1$Cauchy.pval,1))/qchisq(0.5,1)
text(0.5,12, paste("lambda","=",  signif(alpha, digits = 3)) )
dev.off()


###### LOF only
allres<-NULL
for (num in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/meta/eomi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",num,".RData")
load(filename)
sum0<-result[["mask1"]]
sum0$chr<-num
#sum1<-subset(sum0,Group %in% siggenes$Group)
allres<-rbind(allres,sum0)
}
allres$type<-"hclof_noflag"
lof<-allres


####
#### manhattan plot
aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(lof,gencode[,c("geneid","gene_name","start")],by.x="Group",by.y="geneid")

setwd("/medpop/esp2/schoi/migen/meta/figure")

library(qqman)

comb2<-unique(comb1[,c("Group","chr","start","gene_name","Cauchy.anytranscript.pval")])

png("eomi_UKBB_AoU_MIGEN_MGBB_hclof_noflag_Burden_Cauchy_manhattan.png",width = 1000, height = 500,res=110)
manhattan(comb2,chr="chr",bp="start",p="Cauchy.anytranscript.pval",snp="gene_name",col = c("#1b9e77","#d95f02"),genomewideline = -log10(gthreshold),
    suggestiveline = -log10(sthreshold),annotatePval=gthreshold,annotateTop = FALSE)
dev.off()


png("eomi_UKBB_AoU_MIGEN_MGBB_hclof_noflag_Burden_Cauchy_qq.png",width = 1000, height = 1000,res=110)
qq(comb2$Cauchy.anytranscript.pval)
alpha<-median(qchisq(1-comb2$Cauchy.anytranscript.pval,1))/qchisq(0.5,1)
text(0.5,12, paste("lambda","=",  signif(alpha, digits = 3)) )
dev.off()

loflist<-subset(comb2,Cauchy.anytranscript.pval<gthreshold)


######
###### hclof_noflag_missense0.8
allres<-NULL
for (num in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/meta/eomi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",num,".RData")
load(filename)
sum0<-result[["mask2"]]
sum0$chr<-num
#sum1<-subset(sum0,Group %in% siggenes$Group)
allres<-rbind(allres,sum0)
}
allres$type<-"hclof_noflag_missense0.8"
lof_missense<-allres

####
####
aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(lof_missense,gencode[,c("geneid","gene_name","start")],by.x="Group",by.y="geneid")

setwd("/medpop/esp2/schoi/migen/meta/figure")

library(qqman)

comb2<-unique(comb1[,c("Group","chr","start","gene_name","Cauchy.anytranscript.pval")])

png("eomi_UKBB_AoU_MIGEN_MGBB_hclof_noflag_missense0.8_Burden_Cauchy_manhattan.png",width = 1000, height = 500,res=110)
manhattan(comb2,chr="chr",bp="start",p="Cauchy.anytranscript.pval",snp="gene_name",col = c("#1b9e77","#d95f02"),genomewideline = -log10(gthreshold),
    suggestiveline = -log10(sthreshold),annotatePval=gthreshold,annotateTop = FALSE)
dev.off()


png("eomi_UKBB_AoU_MIGEN_MGBB_hclof_noflag_missense0.8_Burden_Cauchy_qq.png",width = 1000, height = 1000,res=110)
qq(comb2$Cauchy.anytranscript.pval)
alpha<-median(qchisq(1-comb2$Cauchy.anytranscript.pval,1))/qchisq(0.5,1)
text(0.5,12, paste("lambda","=",  signif(alpha, digits = 3)) )
dev.off()




####
### mi
source("/medpop/esp2/schoi/software/meta_transcript/src/rarevar_tests.R")
source("/medpop/esp2/schoi/software/meta_transcript/src/transcript_cross_annot_meta.R")

for (num in 1:22){

study_path<-c(paste0("/medpop/esp2/schoi/migen/meta/hclof_noflag/UKBB_AoU_MIGEN_MGBB/mi/mi_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_POPMAX0.001_SPA_chr",num,".RData"),paste0("/medpop/esp2/schoi/migen/meta/hclof_noflag_missense0.8/UKBB_AoU_MIGEN_MGBB/mi/mi_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_missense0.8_POPMAX0.001_SPA_chr",num,".RData"))

result<-transcript_cross_annot_meta(study_path=study_path,test=c("Burden"),
                                   min_study_cmac=20,min_meta_cmac=20,
                                   use.anytranscript=T)

outfile<-paste0("/medpop/esp2/schoi/migen/meta/mi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",num,".RData")
save(result,file=outfile)
}

########
########
R.utils::sourceDirectory("/medpop/esp2/schoi/software/meta_transcript/src")

library(GENESIS)
library(CompQuadFrom)
library(survey)

######
###### manhattan plot
allres<-NULL
for (num in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/meta/mi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",num,".RData")
load(filename)
sum0<-result[["Cauchy.result"]]
sum0$chr<-num
allres<-rbind(allres,sum0)
}
nrow(allres)

####
#### manhattan plot
aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(allres,gencode[,c("geneid","gene_name","start")],by.x="Group",by.y="geneid")
gthreshold<-0.05/nrow(comb1)
sthreshold<-0.1/nrow(comb1)

siggenes<-subset(comb1,Cauchy.pval<gthreshold)
siggenelist<-siggenes$gene_name

setwd("/medpop/esp2/schoi/migen/meta/figure")
library(qqman)

png("mi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_manhattan.png",width = 1000, height = 500,res=110)
manhattan(comb1,chr="chr",bp="start",p="Cauchy.pval",snp="gene_name",col = c("#1b9e77","#d95f02"),genomewideline = -log10(gthreshold),
    suggestiveline = -log10(sthreshold),annotatePval=gthreshold,annotateTop = FALSE)
dev.off()


png("mi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_qq.png",width = 1000, height = 1000,res=110)
qq(comb1$Cauchy.pval)
alpha<-median(qchisq(1-comb1$Cauchy.pval,1))/qchisq(0.5,1)
text(0.5,12, paste("lambda","=",  signif(alpha, digits = 3)) )
dev.off()


###### LOF only
allres<-NULL
for (num in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/meta/mi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",num,".RData")
load(filename)
sum0<-result[["mask1"]]
sum0$chr<-num
#sum1<-subset(sum0,Group %in% siggenes$Group)
allres<-rbind(allres,sum0)
}
allres$type<-"hclof_noflag"
lof<-allres


####
#### manhattan plot
aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(lof,gencode[,c("geneid","gene_name","start")],by.x="Group",by.y="geneid")

setwd("/medpop/esp2/schoi/migen/meta/figure")

library(qqman)

comb2<-unique(comb1[,c("Group","chr","start","gene_name","Cauchy.anytranscript.pval")])

png("mi_UKBB_AoU_MIGEN_MGBB_hclof_noflag_Burden_Cauchy_manhattan.png",width = 1000, height = 500,res=110)
manhattan(comb2,chr="chr",bp="start",p="Cauchy.anytranscript.pval",snp="gene_name",col = c("#1b9e77","#d95f02"),genomewideline = -log10(gthreshold),
    suggestiveline = -log10(sthreshold),annotatePval=gthreshold,annotateTop = FALSE)
dev.off()


png("mi_UKBB_AoU_MIGEN_MGBB_hclof_noflag_Burden_Cauchy_qq.png",width = 1000, height = 1000,res=110)
qq(comb2$Cauchy.anytranscript.pval)
alpha<-median(qchisq(1-comb2$Cauchy.anytranscript.pval,1))/qchisq(0.5,1)
text(0.5,12, paste("lambda","=",  signif(alpha, digits = 3)) )
dev.off()

loflist<-subset(comb2,Cauchy.anytranscript.pval<gthreshold)


######
###### hclof_noflag_missense0.8
allres<-NULL
for (num in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/meta/mi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",num,".RData")
load(filename)
sum0<-result[["mask2"]]
sum0$chr<-num
#sum1<-subset(sum0,Group %in% siggenes$Group)
allres<-rbind(allres,sum0)
}
allres$type<-"hclof_noflag_missense0.8"
lof_missense<-allres

####
####
aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(lof_missense,gencode[,c("geneid","gene_name","start")],by.x="Group",by.y="geneid")

setwd("/medpop/esp2/schoi/migen/meta/figure")

library(qqman)

comb2<-unique(comb1[,c("Group","chr","start","gene_name","Cauchy.anytranscript.pval")])

png("mi_UKBB_AoU_MIGEN_MGBB_hclof_noflag_missense0.8_Burden_Cauchy_manhattan.png",width = 1000, height = 500,res=110)
manhattan(comb2,chr="chr",bp="start",p="Cauchy.anytranscript.pval",snp="gene_name",col = c("#1b9e77","#d95f02"),genomewideline = -log10(gthreshold),
    suggestiveline = -log10(sthreshold),annotatePval=gthreshold,annotateTop = FALSE)
dev.off()


png("mi_UKBB_AoU_MIGEN_MGBB_hclof_noflag_missense0.8_Burden_Cauchy_qq.png",width = 1000, height = 1000,res=110)
qq(comb2$Cauchy.anytranscript.pval)
alpha<-median(qchisq(1-comb2$Cauchy.anytranscript.pval,1))/qchisq(0.5,1)
text(0.5,12, paste("lambda","=",  signif(alpha, digits = 3)) )
dev.off()

scp schoi@login:/medpop/esp2/schoi/migen/meta/figure/mi_UKBB_AoU_MIGEN_MGBB_* /Users/seuchoi/Documents/project/lipid/CAD_rare/





####### eomi first
allres<-NULL
for (num in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/meta/eomi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",num,".RData")
load(filename)
sum0<-result[["Cauchy.result"]]
sum0$chr<-num
allres<-rbind(allres,sum0)
}
nrow(allres)

####
#### manhattan plot
aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(allres,gencode[,c("geneid","gene_name","start")],by.x="Group",by.y="geneid")
gthreshold<-0.05/nrow(comb1)
siggenes<-subset(comb1,Cauchy.pval<gthreshold)
siggenelist<-siggenes$gene_name

aa<-load("/medpop/esp2/schoi/migen/meta/canonical_transcript_info_all.RData")

siggenes<-merge(siggenes,cinfo,by.x="Group",by.y="geneid")
siggenes<-siggenes[order(siggenes$chr,siggenes$start),]

for (gg in 1:nrow(siggenes)){
chr<-siggenes[gg,"chr"]
geneid<-siggenes[gg,"Group"]


###ukbb
ukbblof<-get(load(paste0("/medpop/esp2/schoi/migen/ukbb/result/burden/spa/eomi/hclof_noflag_POPMAX0.001/UKBB_migen_ALL_EOMI_hclof_noflag_POPMAX0.001_chr",chr,"_SPA.RData")))
ukbblofmissense<-get(load(paste0("/medpop/esp2/schoi/migen/ukbb/result/burden/spa/eomi/hclof_noflag_missense0.8_7tools_POPMAX0.001/UKBB_migen_ALL_EOMI_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr",chr,"_SPA.RData")))

#####AoU
aoulof<-get(load(paste0("/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_POPMAX0.001/eomi_hclof_noflag_POPMAX0.001_chr",chr,"_SPA.RData")))
aoulofmissense<-get(load(paste0("/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_missense0.8_POPMAX0.001/eomi_hclof_noflag_missense0.8_POPMAX0.001_chr",chr,"_SPA.RData")))

####MIGEN
migenlof<-get(load(paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_POPMAX0.001/spa/MIGEN_V13.EOMI.hclof_noflag_POPMAX0.001_chr",chr,"_SPA.RData")))
migenlofmissense<-get(load(paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/spa/MIGEN_V13.EOMI.hclof_noflag_missense0.8_POPMAX0.001.chr",chr,"_SPA.RData")))

#####MGBB
mgbblof<-get(load(paste0("/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_POPMAX0.001/MGBB.eomi.hclof_noflag_POPMAX0.001.chr",chr,"_SPA.RData")))
mgbblofmissense<-get(load(paste0("/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_missense0.8_POPMAX0.001/MGBB.eomi.hclof_noflag_missense0.8_POPMAX0.001.chr",chr,"_SPA.RData")))

#### meta 
metares<-get(load(paste0("/medpop/esp2/schoi/migen/meta/eomi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",chr,".RData")))
metalof<-metares$mask1
metalofmissense<-metares$mask2

#######
ukbb1<-subset(ukbblof$results,genename==geneid)
if(nrow(ukbb1)>0){
ukbb1$type<-"hclof_noflag"
}
ukbb2<-subset(ukbblofmissense$results,genename==geneid)
if(nrow(ukbb2)>0){
ukbb2$type<-"hclof_noflag_missense0.8"
}
ukbb3<-rbind(ukbb1,ukbb2)
if(nrow(ukbb3)>0){
ukbb3$cohort<-"UKBB"
}

aou1<-subset(aoulof$results,genename==geneid)
aou1<-aou1[,-16]
if(nrow(aou1)>0){
aou1$type<-"hclof_noflag"
}
aou2<-subset(aoulofmissense$results,genename==geneid)
aou2<-aou1[,-16]
if(nrow(aou2)>0){
aou2$type<-"hclof_noflag_missense0.8"
}
aou3<-rbind(aou1,aou2)
if(nrow(aou3)>0){
aou3$cohort<-"AoU"
}


migen1<-subset(migenlof$results,genename==geneid)
if(nrow(migen1)>0){
migen1$type<-"hclof_noflag"
}
migen2<-subset(migenlofmissense$results,genename==geneid)
if(nrow(migen2)>0){
migen2$type<-"hclof_noflag_missense0.8"
}
migen3<-rbind(migen1,migen2)
if(nrow(migen3)>0){
migen3$cohort<-"MIGEN"
}

mgbb1<-subset(mgbblof$results,genename==geneid)
if(nrow(mgbb1)>0){
mgbb1$type<-"hclof_noflag"
}
mgbb2<-subset(mgbblofmissense$results,genename==geneid)
if(nrow(mgbb2)>0){
mgbb2$type<-"hclof_noflag_missense0.8"
}
mgbb3<-rbind(mgbb1,mgbb2)
if(nrow(mgbb3)>0){
mgbb3$cohort<-"MGBB"
}

meta1<-subset(metalof,Group==geneid)
if(nrow(meta1)>0){
meta1$type<-"hclof_noflag"
}
meta2<-subset(metalofmissense,Group==geneid)
if(nrow(meta2)>0){
meta2$type<-"hclof_noflag_missense0.8"
}
meta3<-rbind(meta1,meta2)
if(nrow(meta3)>0){
meta3$cohort<-"META"
}


comco<-rbind(ukbb3,aou3,migen3,mgbb3)
comco1<-subset(comco,n.alt>=20)
comco1$transcript<-ifelse(comco1$transcript=="all","Pseudo_transcript",comco1$transcript)
comco2<-merge(comco1,siggenes,by.x="genename",by.y="Group")
comco2$nlog10p<--log10(comco2$Burden_SPA.pval)


meta4<-merge(meta3,siggenes,by.x="Group",by.y="Group")
meta4$nlog10p<--log10(meta4$Burden.pval)
meta4$Transcript<-ifelse(meta4$Transcript=="all","Pseudo_transcript",meta4$Transcript)


### used variables
### Cauchy.pval,cohort,type,Transcript,nlog10p
comco2.1<-comco2[,c("gene_name","genename","transcript","Burden_SPA.pval","nlog10p","Cauchy.pval","cohort","type")]
meta4.1<-meta4[,c("gene_name","Group","Transcript","Burden.pval","nlog10p","Cauchy.pval","cohort","type")]
names(meta4.1)<-names(comco2.1)
comco3<-rbind(comco2.1,meta4.1)


#metares<-subset(siggenes,Group==geneid)
maxvalue<--log10(min(comco2$Burden_SPA.pval,metares$Cauchy.pval))
genename<-comco3[1,"gene_name"]

transcripts<-unique(comco2$transcript)
xnum<-length(transcripts)+1

canon_transcript<-comco2[1,"ctranscript"]
anytransrcript<-geneid
extratranscripts<-sort(transcripts[!transcripts %in% c(canon_transcript,anytransrcript)],decreasing=T)

ordered_transcripts<-c(canon_transcript,extratranscripts)


### used variables
### Cauchy.pval,cohort,type,Transcript,nlog10p

png(paste0("/medpop/esp2/schoi/migen/meta/figure/UKBB_AoU_MIGEN_MGBB/",genename,"_cohort_type_transcript_v2_test.png"),width = 2000, height = 1000,res=100)
plot(x=NULL,y=NULL,xlim=c(0.5,xnum+0.5),ylim=c(0,maxvalue+2),xaxt="n",ylab="-log10(pvalue)",xlab="Transcripts",main=genename)
axis(side=1,at=c(1:xnum),labels=c("Cauchy Combined",ordered_transcripts))
mtext(text=c("(Canonical)"),side=1,at=c(2),line=1.7)
abline(h=-log10(0.05),col="gray60",lty=2, lwd=1.5)
#dev.off()


### meta
#metares<-subset(siggenes,Group==geneid)
#points(x=1,y=-log10(comco3$Cauchy.pval[1]),pch=23,bg="purple",col="black",cex=2.5)
points(x=1,y=-log10(comco3$Cauchy.pval[1]),pch=8,col="purple",cex=2.5)
### cohorts
cohorts<-c("UKBB","AoU","MIGEN","MGBB","META")
phcs<-c(21,22,23,24,25)
for (co in 1:length(cohorts)){
co1<-cohorts[co]
### hclof_noflag
comco4<-subset(comco3,cohort==co1 & type=="hclof_noflag")

if(nrow(comco4)>0){

for (ta in 1:length(ordered_transcripts)){
target_transcript<-ordered_transcripts[ta]
tt1<-subset(comco4,transcript==target_transcript)
if(nrow(tt1)>0){
points(x=(ta+1),y=tt1$nlog10p,pch=phcs[co],bg="red3",col="black",cex=2)    
}
}
}

### hclof_noflag_missense0.8
comco4<-subset(comco3,cohort==co1 & type=="hclof_noflag_missense0.8")

if(nrow(comco4)>0){

for (ta in 1:length(ordered_transcripts)){
target_transcript<-ordered_transcripts[ta]
tt1<-subset(comco4,transcript==target_transcript)
if(nrow(tt1)>0){
points(x=(ta+1),y=tt1$nlog10p,pch=phcs[co],bg="blue3",col="black",cex=2)    
}
}
}
}
# legend for cohorts
legend("topleft",legend=c("UKBB","AoU","MIGen","MGBB","Meta"),pch=phcs,col="black",pt.bg="grey30",bty="n",pt.cex=2 )
legend("bottomleft",legend=c("Cauchy Combine","hclof_noflag","hclof_noflag_missense0.9"),pch=c(15,15,15),col=c("purple","red3","blue3"),bty="n",pt.cex=2)

dev.off()

}

scp schoi@login:/medpop/esp2/schoi/migen/meta/figure/UKBB_AoU_MIGEN_MGBB/* /Users/seuchoi/Documents/project/lipid/CAD_rare/

scp schoi@login:/medpop/esp2/schoi/migen/meta/figure/UKBB_AoU_MIGEN_MGBB/NUP98_cohort_type_transcript_v2_test.png /Users/seuchoi/Documents/project/lipid/CAD_rare/


cinfo0<-NULL
for (num in 1:22){

load(paste0("/medpop/esp2/projects/UK_Biobank/WES_450K/deepVariant/annotation_v2/out/grouping/missense/ukb23148.chr",num,".genotype_QCed.merged.pvar.vep105.gz.missense.RData"))

gp1<-subset(group,CANONICAL=="YES")
cinfo1<-data.frame(geneid=gp1$group_id,ctranscript=gp1$TranscriptID)
cinfo1<-unique(cinfo1)
cinfo0<-rbind(cinfo0,cinfo1)
}
cinfo<-cinfo0
save(cinfo,file="/medpop/esp2/schoi/migen/meta/canonical_transcript_info_all.RData")


####
#### 

allres<-NULL
for (num in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/meta/eomi_UKBB_AoU_MIGEN_MGBB_Burden_Cauchy_chr",num,".RData")
load(filename)
sum0<-result[["Cauchy.result"]]
sum0$chr<-num
allres<-rbind(allres,sum0)
}
nrow(allres)

####
#### manhattan plot
aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(allres,gencode[,c("geneid","gene_name","start")],by.x="Group",by.y="geneid")
gthreshold<-0.05/nrow(comb1)
sthreshold<-0.1/nrow(comb1)

siggenes<-subset(comb1,Cauchy.pval<gthreshold)
siggenes2<-siggenes[order(siggenes$chr,siggenes$start),]
siggenelist<-siggenes$gene_name


chr<-siggenes2$chr[1]
geneid<-siggenes2$Group[1]
cohort<-"UKBB"
mask<-"lofmissense"

test<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="UKBB")
test2<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MIGEN")
test3<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MIGEN")


chr<-siggenes2$chr[2]
geneid<-siggenes2$Group[2]

test<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MIGEN")
test2<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MIGEN")

test3<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="UKBB")
test4<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MGBB")

chr<-siggenes2$chr[3]
geneid<-siggenes2$Group[3]

test<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MIGEN")
test2<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MIGEN")

test3<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="UKBB")
subset(test3$variantInfo,variant.id=="5:31406895:A:C")
test4<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MGBB")
subset(test4$variantInfo,variant.id=="5:31406895:A:C")


#### RICTOR
chr<-siggenes2$chr[4]
geneid<-siggenes2$Group[4]

test<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MIGEN")
test2<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MIGEN")

test3<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="UKBB")
subset(test3$variantInfo,variant.id=="5:38950055:TAATAG:T")
test4<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MGBB")
subset(test4$variantInfo,variant.id=="5:38950055:TAATAG:T")


### LOX
chr<-siggenes2$chr[5]
geneid<-siggenes2$Group[5]

test<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MIGEN")
test2<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MIGEN")

test3<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="UKBB")
subset(test3$variantInfo,variant.id=="5:122074046:G:T")
test4<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MGBB")
subset(test4$variantInfo,variant.id=="5:122074046:G:T")

###FGFR1
chr<-siggenes2$chr[6]
geneid<-siggenes2$Group[6]

test<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MIGEN")
test2<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MIGEN")

test3<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="UKBB")
subset(test3$variantInfo,variant.id=="8:38417884:A:C")
test4<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MGBB")
subset(test4$variantInfo,variant.id=="8:38417884:A:C")


###ZBTB43
chr<-siggenes2$chr[7]
geneid<-siggenes2$Group[7]

test<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MIGEN")
test2<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MIGEN")

test3<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="UKBB")
subset(test3$variantInfo,variant.id=="9:126833715:T:G")
test4<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MGBB")
subset(test4$variantInfo,variant.id=="9:126833715:T:G")

#### BASE1
chr<-siggenes2$chr[8]
geneid<-siggenes2$Group[8]

test<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MIGEN")
test2<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MIGEN")

test3<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="UKBB")
test4<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="UKBB")

#### CMAS
chr<-siggenes2$chr[9]
geneid<-siggenes2$Group[9]

test<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MIGEN")
test2<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MIGEN")

test3<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="UKBB")
test4<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="UKBB")


#### SMARCA4
chr<-siggenes2$chr[10]
geneid<-siggenes2$Group[10]

test<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MIGEN")
test2<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MIGEN")

test3<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="UKBB")
test4<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="UKBB")


#### LDLR
chr<-siggenes2$chr[11]
geneid<-siggenes2$Group[11]

test<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MIGEN")
test2<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MIGEN")

test3<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="UKBB")
test4<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="UKBB")


#### EIF3D
chr<-siggenes2$chr[12]
geneid<-siggenes2$Group[12]

test<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="MIGEN")
test2<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="MIGEN")

test3<-info(trait="eomi",chr=chr,geneid=geneid,mask="lofmissense",cohort="UKBB")
test4<-info(trait="eomi",chr=chr,geneid=geneid,mask="lof",cohort="UKBB")







info<-function(trait,chr,geneid,cohort,mask){


if(cohort=="UKBB"){

if(trait=="eomi"){

if(mask=="lof"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/ukbb/result/burden/spa_v2/eomi/hclof_noflag_POPMAX0.001/UKBB_migen_ALL_eomi_hclof_noflag_POPMAX0.001_chr",chr,"_SPA_v2.RData")))
}else if(mask=="lofmissense"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/ukbb/result/burden/spa_v2/eomi/hclof_noflag_missense0.8_7tools_POPMAX0.001/UKBB_migen_ALL_eomi_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr",chr,"_SPA_v2.RData")))
}else{}

}else if(trait=="mi"){
if(mask=="lof"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/ukbb/result/burden/spa_v2/mi/hclof_noflag_POPMAX0.001/UKBB_migen_ALL_mi_hclof_noflag_POPMAX0.001_chr",chr,"_SPA_v2.RData")))
}else if(mask=="lofmissense"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/ukbb/result/burden/spa_v2/mi/hclof_noflag_missense0.8_7tools_POPMAX0.001/UKBB_migen_ALL_MI_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr",chr,"_SPA_v2.RData")))
}else{}

}else if(trait=="cad"){
if(mask=="lof"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/ukbb/result/burden/spa_v2/cad/hclof_noflag_POPMAX0.001/UKBB_migen_ALL_cad_hclof_noflag_POPMAX0.001_chr",chr,"_SPA_v2.RData")))
}else if(mask=="lofmissense"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/ukbb/result/burden/spa_v2/cad/hclof_noflag_missense0.8_7tools_POPMAX0.001/UKBB_migen_ALL_cad_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr",chr,"_SPA_v2.RData")))
}else{}
}
}

if(cohort=="AOU"){

if(trait=="eomi"){

if(mask=="lof"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_POPMAX0.001/eomi_hclof_noflag_POPMAX0.001_chr",chr,"_SPA.RData")))
}else if(mask=="lofmissense"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_missense0.8_POPMAX0.001/eomi_hclof_noflag_missense0.8_POPMAX0.001_chr",chr,"_SPA.RData")))
}else{}

}else if(trait=="mi"){
if(mask=="lof"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_POPMAX0.001/mi_hclof_noflag_POPMAX0.001_chr",chr,"_SPA.RData")))
}else if(mask=="lofmissense"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_missense0.8_POPMAX0.001/mi_hclof_noflag_missense0.8_POPMAX0.001_chr",chr,"_SPA.RData")))
}else{}

}else if(trait=="cad"){
if(mask=="lof"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_POPMAX0.001/cad_hclof_noflag_POPMAX0.001_chr",chr,"_SPA.RData")))
}else if(mask=="lofmissense"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_missense0.8_POPMAX0.001/cad_hclof_noflag_missense0.8_POPMAX0.001_chr",chr,"_SPA.RData")))
}else{}
}
}

if(cohort=="MIGEN"){

if(trait=="eomi"){

if(mask=="lof"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_POPMAX0.001/spa_v2/MIGEN_V13.EOMI.hclof_noflag_POPMAX0.001_chr",chr,"_SPA_v2.RData")))
}else if(mask=="lofmissense"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/spa_v2/MIGEN_V13.EOMI.hclof_noflag_missense0.8_POPMAX0.001.chr",chr,"_SPA_v2.RData")))
}else{}

}else if(trait=="mi"){
if(mask=="lof"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_POPMAX0.001/spa_v2/MIGEN_V13.MI.hclof_noflag_POPMAX0.001_chr",chr,"_SPA_v2.RData")))
}else if(mask=="lofmissense"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/spa_v2/MIGEN_V13.MI.hclof_noflag_missense0.8_POPMAX0.001.chr",chr,"_SPA_v2.RData")))
}else{}

}else if(trait=="cad"){
if(mask=="lof"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_POPMAX0.001/spa_v2/MIGEN_V13.CHD.hclof_noflag_POPMAX0.001_chr",chr,"_SPA_v2.RData")))
}else if(mask=="lofmissense"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/spa_v2/MIGEN_V13.CHD.hclof_noflag_missense0.8_POPMAX0.001.chr",chr,"_SPA_v2.RData")))
}else{}
}
}


if(cohort=="MGBB"){

if(trait=="eomi"){

if(mask=="lof"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_POPMAX0.001/spa_v2/MGBB.eomi.hclof_noflag_POPMAX0.001.chr",chr,"_SPA_v2.RData")))
}else if(mask=="lofmissense"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_missense0.8_POPMAX0.001/spa_v2/MGBB.eomi.hclof_noflag_missense0.8_POPMAX0.001.chr",chr,"_SPA_v2.RData")))
}else{}

}else if(trait=="mi"){
if(mask=="lof"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_POPMAX0.001/spa_v2/MGBB.mi.hclof_noflag_POPMAX0.001.chr",chr,"_SPA_v2.RData")))
}else if(mask=="lofmissense"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_missense0.8_POPMAX0.001/spa_v2/MGBB.mi.hclof_noflag_missense0.8_POPMAX0.001.chr",chr,"_SPA_v2.RData")))
}else{}

}else if(trait=="cad"){
if(mask=="lof"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_POPMAX0.001/spa_v2/MGBB.cad.hclof_noflag_POPMAX0.001.chr",chr,"_SPA_v2.RData")))
}else if(mask=="lofmissense"){
dat<-get(load(paste0("/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_missense0.8_POPMAX0.001/spa_v2/MGBB.cad.hclof_noflag_missense0.8_POPMAX0.001.chr",chr,"_SPA_v2.RData")))
}else{}
}
}



res1<-subset(dat$result,genename==geneid)
res2<-dat$variantInfo[[geneid]]
out<-list()
out[["result"]]<-res1
if(cohort!="AOU"){
out[["variantInfo"]]<-res2
}
return(out)
}





load("MIGen.EOMI.covar.nodup.nullmodel.RData")
aa<-table(nullmod$outcome)
load("MIGen.MI.covar.nodup.nullmodel.RData")
bb<-table(nullmod$outcome)
load("MIGen.CHD.covar.nodup.nullmodel.RData")
cc<-table(nullmod$outcome)

rbind(aa,bb,cc)
       0    1
aa 16187 6187
bb 15946 7292
cc 15946 7886



load("MGBB.eomi.eur.covar.nodup.nullmodel.RData")
aa<-table(nullmod$outcome)
load("MGBB.mi.eur.covar.nodup.nullmodel.RData")
bb<-table(nullmod$outcome)
load("MGBB.cad.eur.covar.nodup.nullmodel.RData")
cc<-table(nullmod$outcome)

rbind(aa,bb,cc)






load("UKBB_migen_ALL.EOMI.covar.txt.gz_nullmodel.RData")
aa<-table(nullmod$fit$outcome)
load("UKBB_migen_ALL.MI.covar.txt.gz_nullmodel.RData")
bb<-table(nullmod$fit$outcome)
load("UKBB_migen_ALL.CAD.covar.txt.gz_nullmodel.RData")
cc<-table(nullmod$fit$outcome)

rbind(aa,bb,cc)



load("eomi_covariates_nullmodel.RData")
aa<-table(nullmod$outcome)
load("MI_covariates_nullmodel.RData")
bb<-table(nullmod$outcome)
load("cad_covariates_nullmodel.RData")
cc<-table(nullmod$outcome)

rbind(aa,bb,cc)




aa<-load("/medpop/esp2/schoi/migen/meta/hclof_noflag_missense0.8/ UKBB_AoU_MIGEN_MGBB/eomi/eomi_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_missense0.8_POPMAX0.001_SPA_chr11.RData")

mgbblof$variantInfo[["ENSG00000128652"]]
2:176169216:C:A
2:176169243:C:A 



####
#### UK biobank
cd /medpop/afib/

use .python-3.8.3
source dxpy/bin/activate

#dx login
dx select "cad_rvas"

cd /medpop/esp2/schoi/migen/ukbb/result/burden/spa_v2


dx download cad_rvas:/schoi/nullmodel/UKBB_migen_ALL.*.covar.txt.gz_nullmodel.RData
dx download cad_rvas:/schoi/association/spa_v2/hclof_noflag_missense0.8/*




aa<-load("/medpop/esp2/schoi/migen/ukbb/result/burden/spa/eomi/hclof_noflag_missense0.8_7tools_POPMAX0.001/UKBB_migen_ALL_EOMI_hclof_noflag_missense0.8_7tools_POPMAX0.001_chr11_SPA.RData")
subset(assoc$result,transcript=="ENST00000397007")

load("/medpop/esp2/schoi/migen/aou/spa/hclof_noflag_missense0.8_POPMAX0.001/eomi_hclof_noflag_missense0.8_POPMAX0.001_chr11_SPA.RData")
subset(assoc$result,transcript=="ENST00000397007")

load("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/spa/MIGEN_V13.EOMI.hclof_noflag_missense0.8_POPMAX0.001.chr11_SPA.RData")
subset(assoc$result,transcript=="ENST00000397007")

in4=/medpop/esp2/schoi/migen/mgbb/result/association/hclof_noflag_missense0.8_POPMAX0.001/MGBB.eomi.hclof_noflag_missense0.8_POPMAX0.001.chr${num}_SPA.RData
out1=/medpop/esp2/schoi/migen/meta/hclof_noflag_missense0.8/eomi_UKBB_AoU_MIGEN_MGBB_Burden_hclof_noflag_missense0.8_POPMAX0.001_SPA_chr${num}.RData




