# ======================================#
#            Input from user            #
# ======================================#

  # pheno_name: starting string of the output
  # pheno_file: name of phenotype file; phenotype file should be in .tsv format. 

  # id_col:      column name for ID in pheno_file
  # outcome_col: column name for the binary outcome (with value 0/1) or continuous outcome in pheno_file
  # covars:      covariates that we want to adjust for the NULL model        
  # covars_test: covariates that we want to test before adjusting in the NULL model (only keep the significant ones)
  # relatedness: sparse GRM  (RData), e.g.,
                 # n x n sparse Matrix of class "dsCMatrix"
                 #        1000000 1000001 1000002  1000003 1000004
                 #1000000 1.00000 0.49000 0.20000  .       0.22222
                 #1000001 0.49000 1.00000 0.20006  0.33333 .
                 #1000002 0.20000 0.20006 1.00000  .       0.44444
                 #1000003 .       0.33333 .        1.0000  .
                 #1000004 0.22222 .       0.44444  .       1.00000

  # unrelated: tsv file includes 1 column with "ID" of unrelated individuals, e.g.,
               # ID
               # 123456
               # 111111
               # 164666

  # chr:       chromosome number (e.g., 12)
  # gene_id:   gene ID (e.g., "ENSG00000155657")
  # gene_name: gene name (e.g., "TTN")  
  # gdsfile:   file name of the genotype data (GDS file)
  # groupfile: file name of the grouping (annotation) information (RData). 
               # The grouping file should be a single dataframe called 'group' that is saved within a .RData file.
               # It should contain the following columns:
                 # chr           : chr  "2" "2" "2" "2" ...
                 # pos           : num  178527023 178527023 178527025 178527025  ...
                 # ref           : chr  "G" "G" "T" "T" ...
                 # alt           : chr  "A" "A" "C" "C" ...
                 # group_id      : chr  "ENSG00000155657" "ENSG00000155657" "ENSG00000155657" "ENSG00000155657" ...
                 # CANONICAL     : chr  "-" "-" "-" "-" ...
                 # TranscriptID  : chr  "ENST00000342175" "ENST00000342992" "ENST00000359218" "ENST00000460472" ...


############################################################
source('/UKBB_200KWES_CVD/GENESIS_adaptation_source.R')
source("/meta_transcript/src/ExtractKernalStatistics_SPA_transcript.R")
#source("/meta_transcript/src/transcript_meta_analysis.R")
source("/meta_transcript/src/CCT.R")

t_aware_analysis <- function(pheno_name, pheno_file, id_col, outcome_col, covars, covars_test, 
                             relatedness=NULL, unrelated=NULL, #out_dir,
                             chr, gene_id=NULL, gene_name=NULL, gdsfile, groupfile,

                             Test_Covar_P_cutoff=0.05, Model_type="binomial"){



      # output names
        outname_null <- paste0(pheno_name, "_nullmodel.RData") # name of the output from NULL model
        #outname_null <- paste0(out_dir, outname_null)

        if(!is.null(gene_id) && !is.null(gene_name)){
	          # (1) For one gene
		          outname_assoc <- paste0(pheno_name, "_chr", chr, "_", gene_name, "_assoc.RData")
  			      outname_final <- paste0(pheno_name, "_chr", chr, "_", gene_name, "_final.RData")
  			
		    }else{ 
		        # (2) For one chromosome
		          outname_assoc <- paste0(pheno_name, "_chr", chr, "_assoc.RData")
		          outname_final <- paste0(pheno_name, "_chr", chr, "_final.RData")
	    	}

      # model type
        if(Model_type == "binomial"){       IV_Rank_Norm <- FALSE
              						                  vc.test      <- "Score.SPA"
        }else if(Model_type == "gaussian"){ IV_Rank_Norm <- TRUE
              						                  vc.test      <- "Score"      
        }else{ stop("Invalid model type") }

        cat("########################################")
        cat(paste0("Running analysis for ", Model_type, " using ", vc.test, ". IV Rank Normalization is: ", IV_Rank_Norm))
        cat("\n Fixed covariates in the null model: ")
        cat(covars, sep = ", ")
       

      # ----------------------- (1)  NULL model ------------------------- 
        fit_nullmodel(phenofile=pheno_file, ID_col=id_col, Outcome=outcome_col, IV_Rank_Norm=IV_Rank_Norm, 
                      Fixed_Covars=covars, Test_Covars=covars_test,
                  	  Test_Covar_P_cutoff=Test_Covar_P_cutoff, Model_type=Model_type,
                  	  relfile=relatedness, unrelfile=unrelated, 
                  	  outfile=outname_null)


      # ----------------------- (2) association test by transcripts ------------
        try(kernell_variance_component_aou(
                gdsfile=gdsfile, groupfile=groupfile, phenfile=pheno_file, ID_col=id_col, 
			          nullfile=outname_null, outfile=outname_assoc,
			          test="ExtractKernelStatistics", vc.test=vc.test,
			          AF.max=0.001, MAC.max=Inf, use.weights=FALSE), silent=F)


	   # ----------------------- (3) aggregating --------------------------
	   	 load(outname_assoc)
       res0     <- assoc$results
       res1     <- subset(res0,n.sample.alt>=20)   
       genes    <- unique(res1$genename)
       pval_col <- if(vc.test == "Score.SPA"){"Burden_SPA.pval"}else{"Burden_Score.pval"}
        
       result1<-NULL
       for(gg in 1:length(genes)){
               gene0   <- genes[gg]
               result2 <- subset(res1,genename==gene0)

               #pvals   <- unique(result2$Burden_SPA.pval)
               pvals   <- unique(result2[[pval_col]])
               newpval <- CCT(pvals)

               #result2 <- data.frame(genename=gene0,pval=newpval)
               #result1 <- rbind(result1,result2)
               result1 <- rbind(result1, data.frame(genename=gene0, pval=newpval))
       }

       cauchy <- result1   
       assoc  <- list("results"=assoc$results, "cauchy"=cauchy)
        
       save(assoc,file=outname_final)
  }

