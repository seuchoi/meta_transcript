# Transcript-Aware Rare Genetic Variant Association Analyses  

Docker file is available at: `us.gcr.io/aou-project-385720/meta_transcript_2.18.0:10282025` (archived `us.gcr.io/aou-project-385720/genesis_meta_2.18.0:Sep172024`)


Pre-installed Rpackages are available

```
wget https://bit.ly/3rmveIc -O rpackages4_1_3_aou.tar.gz
tar -xvf rpackages4_1_3_aou.tar.gz
```

Installing time is typically less than 10 mins.



# Preparation - User Input
+ `pheno_name`   Name of the phenotype you are interested. Will be the starting string of the output file. 
  
+ `pheno_file`   File name of the main data set with all covariates and outcome of interest as columns. Phenotype file should be in .tsv format. 
  
+ `id_col`       Column name for sample ID.
   
+ `outcome_col`  Column name for the binary outcome (with value 0/1) or continuous outcome in pheno_file.
  
+ `covars`       List of covariates that you want to adjust for the NULL model.   
  
+ `covars_test`  Covariates that you want to test before adjusting in the NULL model (only keep the significant ones).
  
+ `Model_type`   "binomial" for binary outcome; "gaussian" for continuous outcome.
  
+ `Test_Covar_P_cutoff` Significance cut-off for `covars_test` in the NULL model (default 0.05).

+ `relatedness` (optional)    n x n sparse Genetic Relationship Matrix (GRM). The file should be saved as a sparse matrix called 'sparseMat' that is saved within a RData file.
  For example, n x n sparse Matrix of class "dsCMatrix" with sample ID as row and column names:

    |         | **1000000** | **1000001** | **1000002** | **1000003** |
    |:-------:|:-------:|:-------:|:-------:|:-------:|  
    | **1000000** | 1.00000 | 0.49000 | 0.20000 | .       |
    | **1000001** | 0.49000 | 1.00000 | 0.20006 | 0.33333 |
    | **1000002** | 0.20000 | 0.20006 | 1.00000 | .       |
    | **1000003** | .       | 0.33333 | .       | 1.0000  |

   This should be used if the user wants to perform separate.residual.variances = A character string specifying the name of a categorical variable in the phenotype file to be used to compute separate residual error variances for heterogeneous groups.
                 
+ `unrelated` (optional) A tsv file includes 1 column "ID" that includes IDs of unrelated individuals.

   | **ID** |
   |:------:|
   | 123456 |
   | 111111 |
   | 164666 |
                
+ `chr`        Chromosome number of the genotype file (e.g., 12).
  
+ `gdsfile`    File name of the genotype data (as GDS file).
  
+ `groupfile` File name of the grouping (annotation) information (as RData). The grouping file should be a single dataframe called 'group' that is saved within a RData file.
               It should contain the following columns with the exact same names:
    + **chr** (string): "2" "2" "2" "2" ...
    + **pos** (numeric): 178527021 178527021 178527023 178527023  ...
    + **ref** (string): "T" "T" "G" "G" ... 
    + **alt** (string): "C" "C" "A" "A" ... 
    + **group_id**   (e.g., Ensembl Gene ID): "ENSG00000155657" "ENSG00000155657" "ENSG00000155657" "ENSG00000155657" ...
    + **CANONICAL**  (string, "-" for non-canonical, "YES" for canonical): chr  "-" "-" "-" "YES" ...
    + **TranscriptID** (string): "ENST00000342175" "ENST00000342992" "ENST00000359218" "ENST00000589042" ...
    + `groupfile` must reflect the specified `gene_id` and `gene_name` (either include annotation information for the whole chromosome or pre-specified genes).

                                                    
+ `gene_id`    (Optional) Ensembl Gene ID (e.g., "ENSG00000155657"). Set to NULL if you want to run the analysis for the whole chromosome;
  Set to a specific gene ID if you only want to run the analysis for a specific gene not whole chromosome.
  Note: If `gene_id` and `gene_name` are available, update `groupfile` to include only annotations for the target gene.
  
+ `gene_name`  (Optional) gene name (e.g., "TTN"). Set to NULL if you want to run the analysis for the whole chromosome;
  Set to a specific gene name if you only want to run the analysis for a specific gene not whole chromosome.
  Note: If `gene_id` and `gene_name` are available, update `groupfile` to include only annotations for the target gene.




# Analysis - Example Code   
```
  source("/meta_transcript/src/t_aware_analysis_function.R")
  pheno_name          <- "Cardiomyopathy"           
  pheno_file          <- "your_phefile.tsv"          
  id_col              <- "ID"                        
  outcome_col         <- "cardiomyopathy_status"     
  covars              <- c("age","sex","PC1","PC2")  
  covars_test         <- c("PC3", "PC4")            
  Model_type          <- "binomial"                 
  Test_Covar_P_cutoff <- 0.05                
  relatedness         <- "MyGRM_matrix.RData"     # or NULL if not available
  unrelated           <- "unrelated_ID.tsv"       # or NULL if not available  
```


## (1) Run analysis on the whole chromosome 12:
```
      chr          <- 12                   
      gdsfile      <- "exome_genotype_QCed_chr12.gds"         
      groupfile    <- "annotation_vep109_chr12.RData"
      gene_id      <- NULL             
      gene_name    <- NULL

      t_aware_analysis(
        pheno_name=pheno_name, pheno_file=pheno_file, id_col=id_col, outcome_col=outcome_col,
        covars=covars, covars_test=covars_test, relatedness=relatedness, unrelated=unrelated,
        chr=chr, gene_id=gene_id, gene_name=gene_name, gdsfile=gdsfile, groupfile=groupfile,
        Test_Covar_P_cutoff=Test_Covar_P_cutoff, Model_type=Model_type
      )
```

## (2) Run analysis on 1 gene (TTN from chromosome 2):
```
      chr          <- 2       
      gdsfile      <- "exome_genotype_QCed_chr2.gds"          
      groupfile    <- "annotation_vep109_TTN.RData"
      gene_id      <- "ENSG00000155657"
      gene_name    <- "TTN"

      t_aware_analysis(
        pheno_name=pheno_name, pheno_file=pheno_file, id_col=id_col, outcome_col=outcome_col,
        covars=covars, covars_test=covars_test, relatedness=relatedness, unrelated=unrelated,
        chr=chr, gene_id=gene_id, gene_name=gene_name, gdsfile=gdsfile, groupfile=groupfile,
        Test_Covar_P_cutoff=Test_Covar_P_cutoff, Model_type=Model_type
      )
```


More example code with simulated data for test runs can be found in `example_data` folder.

Other example code for running transcript-aware rare variant analyses step by step can be found in `src/Analysis_by_steps.R`

Runtime depends on the number of samples and variants.


# Output

The final outcome with all the results is stored in a file named "Pheno_name_chr_final.RData" (e.g., Cardiomyopathy_chr12_final.RData or Cardiomyopathy_chr12_TTN_final.RData). It contains a list of two dataframes:

+ `results`: transcript (including pseudo transcript) specific results. P-values for continuous outcome are in **Burden_Score.pval** column and for binary outcome are in **Burden_SPA.pval**. 

+ `cauchy`: trascript-aware results (1 p-value per gene).
