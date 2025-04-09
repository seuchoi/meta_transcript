#######################################################
#    Step 1. fitting null models without variants     #
#######################################################
# https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git    branch v1.2
source('/UKBB_200KWES_CVD/GENESIS_adaptation_source.R')

# For binary outcome:
  filei       = "yourdatafile.tsv"            # file name of your main data set with all covariates and outcome of interest as columns
  ID_col      = "person_id"                   # column name for sample ID  
  y           = "Disease_status"              # column name for the outcome variable (same for the continuous outcome)
  covariates  = c("age", "sex", "PC1", "PC2") # list of column names for all the covariates that you want to adjust as fixed effect
  test_covars = paste0("PC",3:16)             # list of additional covariates that you want to adjust as fixed effect only if they are significant
  test_cutoff = 0.05                          # p-value cutoff for test_covars
  mtype       = "binomial"                    # "binomial" for binary outcome; "gaussian" for continuous outcome
  m_rel       = "Your_sparse_GRM.RData"       # file name for the Relationship Matrix
  unrelfile   = "unrelated.tsv"               # file name for unrelated individuals
  outputname  = "nullmodel_results.RData"     # file name for the null model output
  
  fit_nullmodel(phenofile=filei, 
                ID_col=ID_col, Outcome=y, IV_Rank_Norm=FALSE,      # IV_Rank_Norm=TRUE for continuous outcome
                Fixed_Covars=covariates, Test_Covars=test_covars,
                Test_Covar_P_cutoff=test_cutoff, 
                Model_type=mtype,
                relfile=m_rel, 
                unrelfile=unrelfile, 
                outfile=outputname) 


  
  
#######################################################
#       Step 2.      Association analysis             #
#######################################################
# From rpackages4_1_3
  library(GENESIS)
  library(CompQuadForm)
  library(survey)
  source("/UKBB_200KWES_CVD/GENESIS_adaptation_source.R")
  source("/meta_transcript/src/ExtractKernalStatistics_SPA_transcript.R")


# Example analysis for chromosome 18  
  gdsfile   = "QCed_exome_genotype_variant_chr18.gds"  # your QCed genotype file 
  groupfile = "Annotation_for_variants.RData"          # annotations for the variants (e.g., from VEP) indicating groupings for variants
  phenfile  = "yourdatafile.tsv"                       # file name of your main data set with all covariates and outcome of interest as columns
  ID_col    = "person_id"                              # column name for sample ID  
  nullfile  = "nullmodel_results.RData"                # name of the output file from fit_nullmodel() function in step 1
  outfile   = "outfile_chr18"                          # Prefix for the name of the output from this step

  kernell_variance_component_aou(gdsfile=gdsfile,
                                     groupfile=groupfile,
                                     phenfile=phenfile,
                                     ID_col=ID_col,
                                     nullfile=nullfile,
                                     outfile=outfile,
                                     test="ExtractKernelStatistics",
                                     vc.test="Score.SPA",   # vc.test="Score" for continuous outcome
                                     AF.max=0.001,
                                     MAC.max=Inf,
                                     use.weights=FALSE)
  

  
  
#######################################################
#          Step 3.   Aggregating results              #
#######################################################  
#source("/meta_transcript/src/transcript_meta_analysis.R") 
source("/meta_transcript/src/CCT.R")  
  
# For chromosome 18  
  load("outfile_chr18.RData")              # output from kernell_variance_component_aou() in step 2
  res0  <- assoc$results
  res1  <- subset(res0,n.sample.alt>=20)   
  genes <- unique(res1$genename)
  
  result1 <- NULL
  for(gg in 1:length(genes)){
    gene0   <- genes[gg]
    result2 <- subset(res1,genename==gene0)
    
    pvals   <- unique(result2$Burden_SPA.pval)
    newpval <- CCT(pvals)
    
    result2 <- data.frame(genename=gene0,pval=newpval)
    result1 <- rbind(result1,result2)
  }
  assoc  <- list("results"=assoc$results, "cauchy"=result1)
  
  save(assoc,file="chr18_analysis_results.RData")
  

  
  
  
