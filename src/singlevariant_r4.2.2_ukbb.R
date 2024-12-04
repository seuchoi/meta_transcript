singleassoc_varinfo_ukbb<-function(gdsfile=gdsfile,varfile=NULL,phenfile=phenfile,ID_col="scanID",nullfile=nullfile,stest=test,outfile=outfile){
  gds <- seqOpen(gdsfile, allow.duplicate=T)
  ##### samples
    phen1<-fread(phenfile,header=T,data.table=F,sep="\t")
        names(phen1)[which(colnames(phen1)==ID_col)]<-"sample.id"
        id_int <- FALSE
        if(class(phen1$sample.id)=='integer'){
                id_int <- TRUE
                class(phen1$sample.id) <- 'character'
        }
        samid0<-phen1$sample.id
    
  
  ######
  ###### QCed variants
  if(!is.null(varfile)){
  varinfo<-variantInfo(gds)
  vardata<-fread(varfile,header=F,data.table=F)
  varinfom<-merge(varinfo,vardata,by.x=c("pos","ref","alt"),by.y=c("V2","V3","V4"))
  varid0<-varinfom$variant.id
  }
  ######
  ###### read gds file
         # Read gds file
        #gds <- seqOpen(gdsfile, allow.duplicate=T)
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

  
  ######
  ###### # construct a SeqVarData object
    seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(combphen2))  
  
  ######
  ###### filter the gdsfile
  if(!is.null(varfile)){
  seqSetFilter(seqData, sample.id=samid0, variant.id=varid0)
  }else{
  seqSetFilter(seqData, sample.id=samid0)  
  }
  ######
  ###### create the iterator
  iterator <- SeqVarBlockIterator(seqData, verbose=TRUE)
  
  ######
  ###### load null model
  nullmod<-get(load(nullfile))
  
  ######
  ###### perfrom test
  assoc <- assocTestSingle(iterator, nullmod, test=stest,verbose=TRUE)
  
  ######
  ###### save files
  save(assoc,file=outfile)
  
  #####
  #####
  sessionInfo()
  quit("no")
  
}
