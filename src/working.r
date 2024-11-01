
setwd("/medpop/esp2/schoi/software/")

source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")
#source("TOPMed_AFib_pipeline/DNANexus/kernell_variance_component_modfied.R")
source("meta_transcript/src/ExtractKernalStatistics_SPA_transcript.R")
#source("TOPMed_AFib_pipeline/DNANexus/ExtractKernelStatistics_error_fixed.R")


gdsfile="/medpop/esp2/skoyama/passing/migen/ws02/60_output/out/MIGEN_V13_sampleQCed_chr22.gds"
groupfile="/medpop/esp2/skoyama/passing/migen/ws01/70_annotation/out/grouping/hclof_noflag_missense0.8_POPMAX0.001/MIGEN_V13_sampleQCed_chr22.pvar.vep105.gz.hclof_noflag_missense0.8_POPMAX0.001.RData"
phenfile="/medpop/esp2/schoi/migen/migen/data/MIGen.covar.nodup.txt"
ID_col="IID"
collapse=FALSE
nullfile="/medpop/esp2/schoi/migen/migen/data/MIGen.covar.nodup.nullmodel.RData"
outfile="/medpop/esp2/schoi/migen/migen/result/association/MIGEN_V13.hclof_noflag_missense0.8_POPMAX0.001.chr22_spa.RData"
AF.max=0.001; MAC.max=Inf; use.weights=FALSE;
vc.test=c("Score.SPA");
test=c("ExtractKernelStatistics")
SAIGEGENEplus_collapse_threshold=10; weight.beta=c(1,1)


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
        # annot<-subset(annot,group_id=="ENSG00000205560")
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








##### testing assocTestAggregate_Sean



gdsobj=iterator; null.model=nullmod; AF.max=AF.max; MAC.max=MAC.max; test=test; vc.test=vc.test; vc.type="regular weighted"; collapse = FALSE; verbose=TRUE; use.weights=F; weight.user=NULL;
SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold; weight.beta=c(1,1)
neig = 200; ntrace = 500;burden.test="Score.SPA"
rho = seq(from = 0, to = 1, by = 0.1);
sparse=TRUE; imputed=FALSE;
male.diploid=TRUE; genome.build=c("hg38");
verbose=TRUE;grp=gr





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



			   assoc <- testVariantSet_Sean(null.model, G=geno, use.weights=use.weights, weights=weight, freq=freq,
                                              test=test, burden.test=burden.test, collapse=collapse, recessive=recessive, recessive.model=recessive.model,
					      var.info=var.info,
                                              vc.test=vc.test, vc.type=vc.type, SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold,
                                              neig = neig, ntrace = ntrace,
                                              rho=rho,grp=grp)



G=geno; use.weights=use.weights; weights=weight; freq=freq;
test=test; burden.test=burden.test; collapse=collapse; #recessive=recessive; #recessive.model=recessive.model;
var.info=var.info;vc.test=vc.test; vc.type=vc.type; SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold;
neig = neig; ntrace = ntrace;rho=rho



null.model=nullmod; G, weights, freq, use.weights=F, var.info,
                            test = c("Burden", "SKAT", "fastSKAT", "SMMAT", "fastSMMAT", "SKATO", "SKAT_SAIGEGENEplus", "ExtractKernelStatistics"),
                            burden.test = c("Score","Score.SPA"), collapse = FALSE, recessive=FALSE, recessive.model = c("strict", "putative"),
                            vc.type = "regular weighted", vc.test=c("Score","Score.SPA"), SAIGEGENEplus_collapse_threshold=10,
                            neig = 200, ntrace = 500,
                            rho = seq(from = 0, to = 1, by = 0.1)





testVariantSet_Sean <- function( nullmod, G, weights, freq, use.weights=F, var.info,
                            test = c("Burden", "SKAT", "fastSKAT", "SMMAT", "fastSMMAT", "SKATO", "SKAT_SAIGEGENEplus", "ExtractKernelStatistics"),
                            burden.test = c("Score","Score.SPA"), collapse = FALSE, recessive=FALSE, recessive.model = c("strict", "putative"),
                            vc.type = "regular weighted", vc.test=c("Score","Score.SPA"), SAIGEGENEplus_collapse_threshold=10,
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
	# SPA not yet supported.... Will implement later
        # SAIGEGENEplus_collapse not yet implemented ... Will work on this later.
        out <- testVariantSet_ExtractKernelStatistics_ScoresAndCovarianceMatrices_Sean(nullmod, G, weights, var.info, neig = Inf, ntrace = Inf, Use.SPA=Use.SPA, freq=freq,
                                                                                       SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold)
    }
    return(out)
}



out <- testVariantSet_ExtractKernelStatistics_ScoresAndCovarianceMatrices_Sean(nullmod, G, weights, var.info, neig = Inf, ntrace = Inf, Use.SPA=Use.SPA, freq=freq,
                                                                                       SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold)
    }
#####
##### 




testVariantSet_ExtractKernelStatistics_ScoresAndCovarianceMatrices_Sean_spa <- function(nullmod, G, weights, var.info, neig = Inf, ntrace = Inf, Use.SPA=F, freq,
                                                                                       SAIGEGENEplus_collapse_threshold=1){

        # Check for use.SPA, which is not yet supported
        if(Use.SPA){
                stop("SPA not yet implemented for ExtractKernelStatistics function. Stopping.")
        }

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
        single_var_out_all <- out
        colnames(V)<-rownames(V)<-rownames(single_var_out_all)<-var.id.name # added

        # Compute SPA adjusted Sum of Variances (SAIGE-GENE, AJHG), we will use this for SKAT test later
        V_tilde <- single_var_out_all$SPA.Score.Variance
        Vsum_tilde <- sum(single_var_out_all$SPA.Score.Variance)
        
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
        out[['single_var_out']] <- single_var_out_all
        out[['covariance_matrix']] <- V
        return(out)
                                                                                       }

# number of variant sites
                  n.site <- length(unique(var.info$variant.id))

                  # number of alternate alleles
                  n.alt <- sum(geno, na.rm=TRUE)

                  # number of samples with observed alternate alleles > 0
                  n.sample.alt <- sum(rowSums(geno, na.rm=TRUE) >= 0.5)







gdsobj=iterator; null.model=nullmod; AF.max=AF.max; MAC.max=MAC.max; test=test; vc.test=vc.test; vc.type="regular weighted"; collapse = FALSE; verbose=TRUE; use.weights=T; weight.user=NULL;
SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold; weight.beta=c(1,1)





AF.max=1, MAC.max=Inf, use.weights=F,
weight.beta=c(1,1), weight.user=NULL,
test=c("Burden", "SKAT", "fastSKAT", "SMMAT", "SKATO", "SKAT_SAIGEGENEplus", "ExtractKernelStatistics"),
                   burden.test=c("Score", "Score.SPA"), collapse=FALSE, recessive=F, recessive.model=c("strict", "putative"),
                   vc.test=c("Score", "Score.SPA"), vc.type="regular weighted", SAIGEGENEplus_collapse_threshold=10,
                   # pval.method=c("davies", "kuonen", "liu"),

                   neig = 200; ntrace = 500;
                   rho = seq(from = 0, to = 1, by = 0.1);
                   sparse=TRUE; imputed=FALSE;
                   male.diploid=TRUE; genome.build=c("hg38");
                   verbose=TRUE





testVariantSet_ExtractKernelStatistics_ScoresAndCovarianceMatrices_Sean <- function(nullmod, G, weights, var.info, neig = Inf, ntrace = Inf, Use.SPA=F, freq,
                                                                                       SAIGEGENEplus_collapse_threshold=1){

        # Check for use.SPA, which is not yet supported
        if(Use.SPA){
                stop("SPA not yet implemented for ExtractKernelStatistics function. Stopping.")
        }

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

	burden_out <- GENESIS:::.testGenoSingleVarScore(burdentilde, G = burden, resid = nullmod$resid, RSS0 = nullmod$RSS0)
        colnames(burden_out) <- paste0("Burden_", colnames(burden_out))
        single_var_out <- GENESIS:::.testGenoSingleVarScore(Gtilde, G = G, resid = nullmod$resid, RSS0 = nullmod$RSS0)
        colnames(V)<-rownames(V)<-rownames(single_var_out)<-var.id.name # added

        out <- list(NULL)
        out[['burden_out']] <- burden_out
        out[['single_var_out']] <- single_var_out
        out[['covariance_matrix']] <- V
        return(out)
}




/medpop/esp2/schoi/migen/migen/script/kernell_variance_component_migen_spa.R


for chrom in {1..22}
do
qsub -N migen_eomi_chr${chrom}_rv_SPA -v chr=${chrom} migen_EOMI_kernel_variance_component_SPA.qsub
done




/medpop/esp2/schoi/migen/migen/script/kernell_variance_component_migen_spa.R


for chrom in {1..22}
do
qsub -N migen_mi_chr${chrom}_rv_SPA -v chr=${chrom} migen_MI_kernel_variance_component_SPA.qsub
done


group

/medpop/esp2/skoyama/passing/migen/ws01/70_annotation/out/grouping/hclof_noflag_POPMAX0.001/MIGEN_V13_sampleQCed_chr${num}.pvar.vep105.gz.hclof_noflag_POPMAX0.001.RData
/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_POPMAX0.001/MIGEN_V13.EOMI.hclof_noflag_POPMAX0.001_chr${num}_SPA.RData



for chrom in {1..22}
do
qsub -N migen_eomi_chr${chrom}_hclof_noflag_SPA -v chr=${chrom} migen_EOMI_kernel_variance_component_hclof_noflag_SPA.qsub
done



/medpop/esp2/skoyama/passing/migen/ws01/70_annotation/out/grouping/hclof_noflag_POPMAX0.001/MIGEN_V13_sampleQCed_chr${num}.pvar.vep105.gz.hclof_noflag_POPMAX0.001.RData
/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_POPMAX0.001/MIGEN_V13.MI.hclof_noflag_POPMAX0.001_chr${num}_SPA.RData



for chrom in {1..22}
do
qsub -N migen_mi_chr${chrom}_hclof_noflag_SPA -v chr=${chrom} migen_MI_kernel_variance_component_hclof_noflag_SPA.qsub
done

##### hclof_noflag_POPMAX0.001
result0<-NULL
for (num in 1:22){

rfile<-paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_POPMAX0.001/MIGEN_V13.MI.hclof_noflag_POPMAX0.001_chr",num,"_SPA.RData")
load(rfile)
result<-assoc$result
result$chr<-num
result0<-rbind(result0,result)
}
outfile<-"/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_POPMAX0.001/MIGEN_V13.MI.hclof_noflag_POPMAX0.001_chrALL_SPA.RData"
result1<-subset(result0,n.sample.alt>=20)
result<-result1
save(result,file=outfile)


result0<-NULL
for (num in 1:22){

rfile<-paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_POPMAX0.001/MIGEN_V13.EOMI.hclof_noflag_POPMAX0.001_chr",num,"_SPA.RData")
load(rfile)
result<-assoc$result
result$chr<-num
result0<-rbind(result0,result)
}
outfile<-"/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_POPMAX0.001/MIGEN_V13.EOMI.hclof_noflag_POPMAX0.001_chrALL_SPA.RData"
result1<-subset(result0,n.sample.alt>=20)
result<-result1
save(result,file=outfile)


#####
#####

##### hclof_noflag_missense0.8_POPMAX0.001
result0<-NULL
for (num in 1:22){

rfile<-paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/MIGEN_V13.MI.hclof_noflag_missense0.8_POPMAX0.001.chr",num,".RData")
load(rfile)
result<-assoc$result
result$chr<-num
result0<-rbind(result0,result)
}
outfile<-"/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/MIGEN_V13.MI.hclof_noflag_missense0.8_POPMAX0.001.chrALL_SPA.RData"
result1<-subset(result0,n.sample.alt>=20)
result<-result1
save(result,file=outfile)


result0<-NULL
for (num in 1:22){

rfile<-paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/MIGEN_V13.EOMI.hclof_noflag_missense0.8_POPMAX0.001.chr",num,"_SPA.RData")

load(rfile)
result<-assoc$result
result$chr<-num
result0<-rbind(result0,result)
}
outfile<-"/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/MIGEN_V13.EOMI.hclof_noflag_missense0.8_POPMAX0.001.chrALL_SPA.RData"
result1<-subset(result0,n.sample.alt>=20)
result<-result1
save(result,file=outfile)

####
####

####
#### CMAS, COLGALT2, GFI1, HOXD3, HECTD4, LOX, NOS3, SMARCA4, ZBTB43

####
####
R.utils::sourceDirectory("/medpop/esp2/schoi/software/meta_transcript/src")

library(GENESIS)
library(CompQuadFrom)
library(survey)

allres<-NULL
for (num in 1:22){

filename<-paste0("/medpop/esp2/schoi/migen/meta/UKBB_AoU_MIGEN_Burden_Cauchy_chr",num,".RData")
load(filename)
sum0<-result[["Cauchy.result"]]
sum0$chr<-num
allres<-rbind(allres,sum0)
}
aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(allres,gencode[,c("geneid","gene_name","start")],by.x="Group",by.y="geneid")
gthreshold<-0.05/nrow(comb1)
sthreshold<-0.1/nrow(comb1)

siggenes<-subset(comb1,Cauchy.pval<gthreshold)
siggenelist<-siggenes$gene_name



#### EOMI for MIGEN hclof_noflag_POPMAX0.001
#### manhattan plot
outfile<-"/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_POPMAX0.001/MIGEN_V13.EOMI.hclof_noflag_POPMAX0.001_chrALL_SPA.RData"
load(outfile)

aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(result,gencode[,c("geneid","gene_name","start")],by.x="genename",by.y="geneid")
gthreshold<-0.05/nrow(comb1)
sthreshold<-0.1/nrow(comb1)

setwd("/medpop/esp2/schoi/migen/meta/figure")

library(qqman)

png("MIGEN_Burden_EOMI_hclof_noflag_POPMAX0.001_Transcripts_manhattan.png",width = 1000, height = 500,res=110)
manhattan(comb1,chr="chr",bp="start",p="Burden_SPA.pval",snp="gene_name",col = c("#1b9e77","#d95f02"),genomewideline = -log10(gthreshold),
    suggestiveline = -log10(sthreshold),annotatePval=gthreshold,annotateTop = FALSE)
dev.off()


png("MIGEN_Burden_EOMI_hclof_noflag_POPMAX0.001_Transcripts_qq.png",width = 1000, height = 1000,res=110)
qq(comb1$Burden_SPA.pval)
alpha<-median(qchisq(1-comb1$Burden_SPA.pval,1))/qchisq(0.5,1)
text(0.5,12, paste("lambda","=",  signif(alpha, digits = 3)) )
dev.off()

#####
##### EOMI for MIGEN hclof_noflag_missense0.8_POPMAX0.001
outfile<-"/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/MIGEN_V13.EOMI.hclof_noflag_missense0.8_POPMAX0.001.chrALL_SPA.RData"
load(outfile)

aa<-load("/medpop/esp2/schoi/migen/annotation/data/gencode.v42.Rdata")
comb1<-merge(result,gencode[,c("geneid","gene_name","start")],by.x="genename",by.y="geneid")
gthreshold<-0.05/nrow(comb1)
sthreshold<-0.1/nrow(comb1)

setwd("/medpop/esp2/schoi/migen/meta/figure")

library(qqman)

png("MIGEN_Burden_EOMI_hclof_noflag_missense0.8_POPMAX0.001_Transcripts_manhattan.png",width = 1000, height = 500,res=110)
manhattan(comb1,chr="chr",bp="start",p="Burden_SPA.pval",snp="gene_name",col = c("#1b9e77","#d95f02"),genomewideline = -log10(gthreshold),
    suggestiveline = -log10(sthreshold),annotatePval=gthreshold,annotateTop = FALSE)
dev.off()


png("MIGEN_Burden_EOMI_hhclof_noflag_missense0.8_POPMAX0.001_Transcripts_qq.png",width = 1000, height = 1000,res=110)
qq(comb1$Burden_SPA.pval)
alpha<-median(qchisq(1-comb1$Burden_SPA.pval,1))/qchisq(0.5,1)
text(0.5,12, paste("lambda","=",  signif(alpha, digits = 3)) )
dev.off()


#####
#####
g1<-subset(gencode,gene_name=="COLGALT2")
chr<-g1$seqnames[1]
geneid<-g1$geneid
migenlofmissense<-get(load(paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/MIGEN_V13.hclof_noflag_missense0.8_POPMAX0.001.chr",chr,".RData")))
res0<-migenlofmissense$variantInfo[[geneid]]


grouping_path<-paste0("/medpop/esp2/skoyama/passing/migen/ws01/70_annotation/out/grouping/hclof_noflag_missense0.8_POPMAX0.001/MIGEN_V13_sampleQCed_chr",chr,".pvar.vep105.gz.hclof_noflag_missense0.8_POPMAX0.001.RData")
load(grouping_path)
target.variants<-subset(group,TranscriptID=="ENST00000421419")

#subset(ukbblofmissense$by.transcript,Transcript==

load(paste0("/medpop/esp2/schoi/migen/migen/result/association/hclof_noflag_missense0.8_POPMAX0.001/MIGEN_V13.hclof_noflag_missense0.8_POPMAX0.001.chr",chr,".RData"))
subset(group,chr==1 & pos==183945462 & ref=="G" & alt=="C")

targetgene<-assoc$variantInfo[geneid]


            variant.id chr       pos   ref alt allele.index n.obs         freq
1      1:183940696:T:C   1 183940696     T   C            1 21426 2.333613e-05
2      1:183944195:C:T   1 183944195     C   T            1 21398 2.336667e-05
3      1:183944314:G:A   1 183944314     G   A            1 21419 2.334376e-05
4      1:183945462:G:C   1 183945462     G   C            1 21363 8.659832e-04
5      1:183945509:G:A   1 183945509     G   A            1 21423 2.333940e-05
6  1:183945544:AGCTG:A   1 183945544 AGCTG   A            1 21422 2.334049e-05
7     1:183945562:G:GC   1 183945562     G  GC            1 21421 2.334158e-05
8      1:183963980:G:T   1 183963980     G   T            1 21309 1.407856e-04
9      1:183975132:G:A   1 183975132     G   A            1 21288 2.348741e-05
10     1:183975147:G:A   1 183975147     G   A            1 21231 2.355047e-05
   MAC weight       Score  Score.SE Score.Stat   Score.pval       Est    Est.SE
1    1      1  0.42475850 0.4942696  0.8593660 3.901386e-01  1.738658 2.0231872
2    1      1 -0.06768964 0.2518201 -0.2688016 7.880824e-01 -1.067435 3.9710891
3    1      1 -0.65383673 0.4757155 -1.3744280 1.693089e-01 -2.889180 2.1020966
4   37      1 19.70389660 1.0747432 18.3335860 4.464755e-75 17.058574 0.9304548
5    1      1  0.50128263 0.4999689  1.0026275 3.160406e-01  2.005380 2.0001242
6    1      1  0.47462449 0.4993224  0.9505372 3.418394e-01  1.903654 2.0027142
7    1      1 -0.25647787 0.4366125 -0.5874268 5.569171e-01 -1.345419 2.2903606
8    6      1 -1.33843538 0.8410957 -1.5912998 1.115421e-01 -1.891937 1.1889254
9    1      1 -0.50579412 0.4998779 -1.0118353 3.116168e-01 -2.024165 2.0004885
10   1      1 -0.03185914 0.1708695 -0.1864531 8.520895e-01 -1.091202 5.8524204
            PVE
1  3.326316e-05
2  3.254400e-06
3  8.508465e-05
4  1.513917e-02
5  4.527792e-05
6  4.069542e-05
7  1.554228e-05
8  1.140542e-04
9  4.611338e-05
10 1.565836e-06






###ukbb
ukbblof<-get(load(paste0("/medpop/esp2/schoi/migen/ukbb/result/burden/hclof_noflag/CAD_HARD_Burden_hclof_noflag_POPMAX0.001_chr",chr,".RData")))
ukbblofmissense<-get(load(paste0("/medpop/esp2/schoi/migen/ukbb/result/burden/hclof_noflag_missense0.8/CAD_HARD_Burden_hclof_noflag_missense0.8_POPMAX0.001_chr",chr,".RData")))

#####AoU
aoulof<-get(load(paste0("/medpop/esp2/schoi/migen/aou/hclof_noflag/AoU_Burden_hclof_noflag_POPMAX0.001_chr",chr,".RData")))
aoulofmissense<-get(load(paste0("/medpop/esp2/schoi/migen/aou/hclof_noflag_missense0.8/AoU_Burden_hclof_noflag_missense0.8_POPMAX0.001_chr",chr,".RData")))

####MIGEN
migenlof<-get(load(paste0("/medpop/esp2/schoi/migen/migen/result/association/burden/hclof_noflag/MIGEN_V13_Burden_hclof_noflag_POPMAX0.001_chr",chr,".RData")))
migenlofmissense<-get(load(paste0("/medpop/esp2/schoi/migen/migen/result/association/burden/hclof_noflag_missense0.8/MIGEN_V13_Burden_hclof_noflag_missense0.8_POPMAX0.001_chr",chr,".RData")))






#####
##### mgbb null model
### UKBB_200KWES_CVD
setwd("/medpop/esp2/schoi/migen/migen/data")
source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")

### 
kinship<-"/medpop/esp2/schoi/migen/migen/data/sample_MGBB_53K_sparse_kinship.RData"
duplicated<-"/medpop/esp2/schoi/migen/mgbb/data/sample_MGBB_53K_sampleQCed.chrall_hg_pruned.segment.all.bed.con"
unrelated<-"/medpop/esp2/schoi/migen/mgbb/data/sample_MGBB_53K_unrelated.tsv"
gdsfile<-"/medpop/esp2/skoyama/passing/mgb_53k_exome_qc/60_output/out/MGBB_53K_sampleQCed.chr22.gds"
phenfile<-"/medpop/esp2/schoi/migen/mgbb/data/MGBB_CAD_inpatient.tsv"

library(data.table)
library(SeqArray)
gds<-seqOpen(gdsfile)
sample.id<-seqGetData(gds,"sample.id")
phen<-fread(phenfile,header=T,data.table=F,sep="\t")
dup<-fread(duplicated,header=T,data.table=F,sep="\t")
phen2<-phen[!phen$Subject_Id %in% dup$ID1,]
phen3<-phen2[phen2$Subject_Id %in% sample.id,]
names(phen3)[1:2]<-c("IID","CAD")
write.table(phen3,"/medpop/esp2/schoi/migen/mgbb/data/MIGen.nocov.nodup.txt",col.names=T,row.names=F,quote=F,sep="\t")

phenfile="/medpop/esp2/schoi/migen/migen/data/MIGen.covar.nodup.txt"


fit_nullmodel(phenofile=phenfile, ID_col="IID", Outcome="CAD", IV_Rank_Norm=FALSE, 
			  Fixed_Covars=c("AGE_BASELINE","SEX"), Test_Covars=paste0("PC",1:10), Test_Covar_P_cutoff=0.05,
			  Model_type=c("binomial"), relfile=kinship, separate.residual.variances=NULL, unrelfile=unrelated, outfile="MIGen.covar.nodup.nullmodel.RData")







.libPaths(c("rpackages4_1_3",.libPaths()))
library(GENESIS)
library(CompQuadForm)
library(survey)
R.utils::sourceDirectory("/home/jupyter/workspaces/cadrvas/meta_transcript/src")

source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")
source("meta_transcript/src/ExtractKernalStatistics_SPA_transcript.R")
num=20


gdsfile<-paste0("/home/jupyter/workspaces/cadrvas/data/gds/exome_genotype_variant_sample_QCed_chr",num,".gds")
groupfile<-paste0("/home/jupyter/workspaces/cadrvas/annotation/AoU_250K_exome_annotation_VEP105_chr",num,".vcf.gz.hclof_noflag_missense0.8_7tools_POPMAX0.001.RData")
phenfile<-"/home/jupyter/workspaces/cadrvas/data/phenotype/CAD_189k_covariates.tsv"
ID_col<-"ID"
nullfile<-"/home/jupyter/workspaces/cadrvas/data/nullmodel/CAD_189k_covariates_nullmodel.RData"
outfile<-paste0("/home/jupyter/workspaces/cadrvas/result/burden/spa/hclof_noflag_missense0.8_POPMAX0.001/CAD_185k_chr",num,"_SPA_Test.RData")

kernell_variance_component_aou(gdsfile=gdsfile,groupfile=groupfile,phenfile=phenfile,ID_col=ID_col,nullfile=nullfile,outfile=outfile, test="ExtractKernelStatistics", vc.test="Score.SPA", AF.max=0.001, MAC.max=Inf,use.weights=FALSE)



setwd("/rprojectnb/adsp-charge/seuchoi/software/")
library(GENESIS)
library(CompQuadForm)
library(survey)
R.utils::sourceDirectory("meta_transcript/src/")
source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")
source("meta_transcript/src/ExtractKernalStatistics_SPA_transcript.R")


num=20
gdsfile<-paste0("/restricted/projectnb/adsp-charge/data/wgs/NG00067_v10_Jul_2023/v0/updatedids_filtered_gcad.qc.compact_filtered.r4.wgs.36361.GATK.2023.06.06.biallelic.genotypes.chr",num,".ALL.gds")
groupfile<-paste0("/restricted/projectnb/adsp-charge/seuchoi/annotation/data/vep_data/grouping/hclof_noflag_POPMAX0.001/updatedids_filtered_gcad.qc.compact_filtered.r4.wgs.36361.GATK.2023.06.06.biallelic.genotypes.chr",num,".ALL.ID.txt.exonbodyfile.annotated.vcf.gz.hclof_noflag_POPMAX0.001.RData")
phenfile<-"/rprojectnb/adsp-charge/data/wgs/36k/phenotype/36k_phenotypes_02202024.tsv"
ID_col<-"SampleID"
nullfile<-"/restricted/projectnb/adsp-charge/seuchoi/36k/result/nullmodel_adsp_36k_r4.0.RData"
outfile<-paste0("/restricted/projectnb/adsp-charge/seuchoi/36k/result/association/adsp_36k_chr",num,"_hclof_noflag_POPMAX0.001.RData")
AF.max=0.001; MAC.max=Inf; use.weights=FALSE;
vc.test=c("Score.SPA");
test=c("ExtractKernelStatistics")
SAIGEGENEplus_collapse_threshold=10; weight.beta=c(1,1)







gdsfile<-paste0("/home/jupyter/workspaces/cadrvas/data/gds/exome_genotype_variant_sample_QCed_chr",num,".gds")
groupfile<-paste0("/home/jupyter/workspaces/cadrvas/annotation/AoU_250K_exome_annotation_VEP105_chr",num,".vcf.gz.hclof_noflag_missense0.8_7tools_POPMAX0.001.RData")
phenfile<-"/home/jupyter/workspaces/cadrvas/data/phenotype/CAD_189k_covariates.tsv"
ID_col<-"ID"
nullfile<-"/home/jupyter/workspaces/cadrvas/data/nullmodel/CAD_189k_covariates_nullmodel.RData"
outfile<-paste0("/home/jupyter/workspaces/cadrvas/result/burden/spa/hclof_noflag_missense0.8_POPMAX0.001/CAD_185k_hclof_noflag_missense0.8_POPMAX0.001_chr",num,"_SPA.RData")
AF.max=0.001; MAC.max=Inf; use.weights=FALSE;
vc.test=c("Score.SPA");
test=c("ExtractKernelStatistics")
SAIGEGENEplus_collapse_threshold=10; weight.beta=c(1,1)


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
        annot <- as.data.frame(annot)
        #annot<-subset(annot,group_id=="ENSG00000130939")
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
                                                 SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold, weight.beta=c(1,1),gr=gr[1:10])
        }else{
                assoc <- assocTestAggregate_Sean(iterator, nullmod, AF.max=AF.max, MAC.max=MAC.max, test=test, vc.test=vc.test, vc.type=vc.type, collapse = FALSE, verbose=TRUE, use.weight=F,
                                                 SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold, weight.beta=weight.beta,gr=gr)
        }












        transcript = ENST00000470736
        genename = ENSG00000130939

        annot<-subset(annot,group_id=="ENSG00000130939")







for(num in 1:22){

groupfile<-paste0("/home/jupyter/workspaces/cadrvas/annotation/AoU_250K_exome_annotation_VEP105_chr",num,".vcf.gz.hclof_noflag_missense0.8_7tools_POPMAX0.001.RData")

load(groupfile)
print(num)
print(names(group))
}






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
    pval <- 1-pcauchy(cct.stat)r
  }
  return(pval)
}



load("PheCode_513.32_count1_hclof_noflag_POPMAX0.001_chr22.RData")
res0<-assoc$results
res1<-subset(res0,n.sample.alt>=20)

genes<-unique(res1$genename)

result1<-NULL
for (gg in 1:length(genes)){
gene0<-genes[gg]
result2<-subset(res1,genename==gene0)

pvals<-unique(result2$Burden_SPA.pval)
newpval<-CCT(pvals)

result2<-data.frame(genename=gene0,pval=newpval)
result1<-rbind(result1,result2)
}





CCT

freqz<-freq$freq



meanImpute_Sean<-function(geno, freqz) {
        nrowz <- nrow(geno)
        ncolz <- ncol(geno)
        try(dimz <- nrowz*ncolz, silent=T)
        if(!is.na(dimz)){
                if(dimz < 1500000000){
                        miss.idx <- Matrix::which(is.na(geno))
                        miss.var.idx <- ceiling(miss.idx/nrow(geno))
                        imputed <- 2*freqz[miss.var.idx]
                        geno[miss.idx] <- 2*freqz[miss.var.idx]
                }else{
                        for(jk in c(1:ncol(geno))){
                                #cat('Busy with', jk, 'out of', ncol(geno), '...\n')
                                miss.rowz <- Matrix::which(is.na(geno[,jk]))
                                if(length(miss.rowz)>0){
                                        imputed <- 2*freqz[jk]
                                        geno[miss.rowz,jk] <- imputed
                                }
                        }
                }
        }else{
                        for(jk in c(1:ncol(geno))){
                                #cat('Busy with', jk, 'out of', ncol(geno), '...\n')
                                miss.rowz <- Matrix::which(is.na(geno[,jk]))
                                if(length(miss.rowz)>0){
                                        imputed <- 2*freqz[jk]
                                        geno[miss.rowz,jk] <- imputed
                                }
                        }
        }
        geno
}