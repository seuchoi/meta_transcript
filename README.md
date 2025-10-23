# Transcript-Aware Rare Genetic Variant Association Analyses  
Example code for running transcript-aware rare variant analyses is in `Analysis_by_steps.R`


Docker file is available at: `us.gcr.io/aou-project-385720/genesis_meta_2.18.0:Sep172024`


Pre-installed Rpackages are available

```
wget https://bit.ly/3rmveIc -O rpackages4_1_3_aou.tar.gz
tar -xvf rpackages4_1_3_aou.tar.gz
```





# Preparation - User Input

```
  pheno_name   <- "Cardiomyopathy"                         # Name of the phenotype we are interested. Will be the starting string of the output file.
  pheno_file   <- "your_phefile.tsv"                       # file name of the main data set with all covariates and outcome of interest as columns. Phenotype file should be in .tsv format. 
  id_col       <- "ID"                                     # column name for sample ID  
  outcome_col  <- "cardiomyopathy_status"                  # column name for the binary outcome (with value 0/1) or continuous outcome in pheno_file
  covars       <- c("ageatdna","genetic_sex","PC1","PC2")  # covariates that we want to adjust for the NULL model.        
  covars_test  <- c("PC3", "PC4")                          # covariates that we want to test before adjusting in the NULL model (only keep the significant ones).
  Model_type   <- "binomial"                               # "binomial" for binary outcome; "gaussian" for continuous outcome.
  Test_Covar_P_cutoff <- 0.05                              # significance cut-off point for covars_test in the NULL model.

  relatedness  <-  "MyGRM_matrix.RData"                    # sparse GRM  (RData), e.g.,
                                                               # n x n sparse Matrix of class "dsCMatrix"
                                                               #        1000000 1000001 1000002  1000003 1000004
                                                               #1000000 1.00000 0.49000 0.20000  .       0.22222
                                                               #1000001 0.49000 1.00000 0.20006  0.33333 .
                                                               #1000002 0.20000 0.20006 1.00000  .       0.44444
                                                               #1000003 .       0.33333 .        1.0000  .
                                                               #1000004 0.22222 .       0.44444  .       1.00000
  unrelated    <- "unrelated_ID.tsv"                       # tsv file includes 1 column with "ID" of unrelated individuals, e.g.,
                                                               # ID
                                                               # 123456
                                                               # 111111
                                                               # 164666


  chr          <- 12                                 # chromosome number (e.g., 12)
  gene_id      <- NULL                               # (Optional) Ensembl Gene ID (e.g., "ENSG00000155657"). Set to a specific gene ID if we only want to run the analysis for a specific gene not whole chromosome.
  gene_name    <- NULL                               # (Optional) gene name (e.g., "TTN"). Set to a specific gene name if we only want to run the analysis for a specific gene not whole chromosome.
  
  gdsfile      <- "exome_genotype_QCed_chr12.gds"    # file name of the genotype data (GDS file)
  groupfile    <- "annotation_vep109_chr12.RData"    # file name of the grouping (annotation) information (RData). 
                                                     # The grouping file should be a single dataframe called 'group' that is saved within a .RData file.
                                                     # It should contain the following columns:
                                                       # chr           : chr  "2" "2" "2" "2" ...
                                                       # pos           : num  178527023 178527023 178527025 178527025  ...
                                                       # ref           : chr  "G" "G" "T" "T" ...
                                                       # alt           : chr  "A" "A" "C" "C" ...
                                                       # group_id      : chr  "ENSG00000155657" "ENSG00000155657" "ENSG00000155657" "ENSG00000155657" ...
                                                       # CANONICAL     : chr  "-" "-" "-" "-" ...
                                                       # TranscriptID  : chr  "ENST00000342175" "ENST00000342992" "ENST00000359218" "ENST00000460472" ...

                                                     # group_id will be the Ensembl Gene ID. 
 
```


# Run the analysis

```
  source("src/t_aware_analysis_function.R")

# (1) Analysis on the whole chromosome 12:
      gene_id      <- NULL             
      gene_name    <- NULL              
      gdsfile      <- "exome_genotype_QCed_chr12.gds"         
      groupfile    <- "annotation_vep109_chr12.RData"

      t_aware_analysis(
        pheno_name=pheno_name, pheno_file=pheno_file, id_col=id_col, outcome_col=outcome_col,
        covars=covars, covars_test=covars_test, relatedness=relatedness, unrelated=unrelated,
        chr=chr, gene_id=gene_id, gene_name=gene_name, gdsfile=gdsfile, groupfile=groupfile,
        Test_Covar_P_cutoff=Test_Covar_P_cutoff, Model_type=Model_type
      )

# (2) Analysis on 1 gene (TTN from chromosome 2):
      gene_id      <- "ENSG00000155657"
      gene_name    <- "TTN"             
      gdsfile      <- "exome_genotype_QCed_chr2.gds"          
      groupfile    <- "annotation_vep109_chr2.RData"

      t_aware_analysis(
        pheno_name=pheno_name, pheno_file=pheno_file, id_col=id_col, outcome_col=outcome_col,
        covars=covars, covars_test=covars_test, relatedness=relatedness, unrelated=unrelated,
        chr=chr, gene_id=gene_id, gene_name=gene_name, gdsfile=gdsfile, groupfile=groupfile,
        Test_Covar_P_cutoff=Test_Covar_P_cutoff, Model_type=Model_type
      )
```
