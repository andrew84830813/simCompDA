


simAddLogNormal = function(n1 = 30, dms_ = 75,seed= 08272008,varParm = 2,meanVec = 0){
  set.seed(seed)
  mu_ = rep(meanVec,dms_)
  sigma_ = diag(dms_)
  U = matrix(runif(dms_*dms_,0,varParm),nrow = dms_)
  U_ = sigma_ + U
  U_ = nearPD(U_)$mat
  eg = min(eigen(as.matrix(sigma_))$values)
  sig = min(eigen(as.matrix(U_))$values)
  dd = min(eg,sig)+0.05
  sig1 = nearPD( U_ + sigma_ + dd*diag(dms_) )$mat
  s1 = sampleAddLogisticNormal(n = n1,dims_ = dms_,mu = mu_ ,sigma = sig1,sigmaScale = 1,sampleName = "S1")
  s1
}



sampleAddLogisticT <-
  function(n, df_,dms,sigmaScale=1,sampleName = "S1",seed= 08272008){
    set.seed(seed)
    sigma = diag(dms)*sigmaScale
    y1 = mvtnorm::rmvt(n,sigma ,df = df_,type = "Kshirsagar")
    y1 = alrInv(y1)
    s1  = data.frame(Status = sampleName,y1)
    return(s1)
  }


simDirDistr <-
  function( n1 = 30 ,dms_ = 75,seed= 08272008,scale = 1){
    
    ## define and scale alpa
    a = rep(1,dms_)*scale
    set.seed(seed)
    s1 = sampleDirichlet(n1 = n1,dims = dms_,sampleName = "S1",a1 = a)$Sample
    
    return(s1)
    
  }





randomKNN_Graph = function(samples,k_){
  k_nn= 1:k_
  adj = diag((samples))
  for (r in 1:nrow(adj)){
    tt = sample(1:ncol(adj))
    tt= ifelse(tt%in%k_nn,1,0)
    adj[r,]=tt
  }
  #Make Matric Symetric
  adj =  knnADJtoSYM(adj)
  
  #Create Graph
  w.graph = graph.adjacency(adj,mode="undirected",weighted = TRUE)
  w.graph.simplified = igraph::simplify(w.graph, remove.loops = TRUE,
                                        edge.attr.comb = igraph_opt("edge.attr.comb"))
  return( w.graph.simplified)
}

s4 = function( n1 = 30 , n2 = 30, dms_ = 75,seed,shift = 2){
  
  
  set.seed(seed)
  set.seed(seed)
  mu_ = rep(0,dms_)
  t25 = round(dms_*.25)
  mu2 = rep(0,dms_)
  mu2[1:t25] = shift
  sigma_ = diag(dms_)
  
  
  U = matrix(runif(dms_*dms_,0,1),nrow = dms_)
  U_ = sigma_ + U
  U_ = nearPD(U_)$mat
  eg = min(eigen(as.matrix(sigma_))$values)
  sig = min(eigen(as.matrix(U_))$values)
  dd = min(eg,sig)+0.05
  sig1 = nearPD( sigma_ + dd*diag(dms_) )$mat
  sig2 = nearPD( sigma_+ U + dd*diag(dms_) )$mat
  
  s1 = sampleAddLogisticNormal(n = n1,dims_ = dms_,mu = mu_ ,sigma = sig1,sigmaScale = 1,sampleName = "S1")
  s2 = sampleAddLogisticNormal(n = n2,dims_ = dms_,mu = mu2,sigma = sig1,sigmaScale = 1,sampleName = "S2")
  df = rbind(s1,s2)
}



getData2 = function(fname,returnCounts = T,binaryConditions=NULL, positiveClass,site = "stool")  {
  path_ = paste(fname,".metaphlan_bugs_list.",site,sep = "")
  loman <- curatedMetagenomicData(path_, dryrun = FALSE)
  loman.eset <- loman[[1]]
  loman.pseq = ExpressionSet2phyloseq( loman.eset )
  HanniganGD.df = otu_table( loman.pseq )@.Data
  HanniganGD.taxa = tax_table( loman.pseq )@.Data
  HanniganGD.df = cbind.data.frame(HanniganGD.taxa,HanniganGD.df)
  HanniganGD.metaData = pData( loman.eset ) 
  data.tree <- ExpressionSet2phyloseq( loman.eset, phylogenetictree = TRUE)
  subID = rownames(HanniganGD.metaData)
  HanniganGD.metaData = HanniganGD.metaData %>% 
    mutate(subjectID = subID) %>% 
    filter(!is.na(study_condition)) 
  
  cnames = colnames(HanniganGD.df)[9:ncol(HanniganGD.df)]
  
  if(returnCounts==T){
    #transform to counts
    for(i in 1:length(cnames)){
      k = which(HanniganGD.metaData$subjectID==cnames[i])
      nreads = HanniganGD.metaData$number_reads[k]
      col = which(colnames(HanniganGD.df)==cnames[i])
      HanniganGD.df[,col] = round(nreads*(HanniganGD.df[,col]/100))
    }
  }
  
  
  HanniganGD.df = HanniganGD.df %>% 
    filter(!is.na(Species)) %>% 
    filter(is.na(Strain)) 
  cn = colnames(HanniganGD.df)[9:ncol(HanniganGD.df)]
  HanniganGD.df  = data.frame(Species = HanniganGD.df$Species,HanniganGD.df[,9:ncol(HanniganGD.df)])
  sp  =as.character(HanniganGD.df$Species)
  HanniganGD.df = data.frame(t(HanniganGD.df[,-1]))
  colnames(HanniganGD.df)=sp
  HanniganGD.df = data.frame(Sample = rownames(HanniganGD.df),HanniganGD.df)
  HanniganGD.df$Sample = as.character(cn)
  md = data.frame(Sample = HanniganGD.metaData$subjectID,Status = HanniganGD.metaData$study_condition)
  md$Sample = as.character(md$Sample)
  HanniganGD.df = left_join(md,HanniganGD.df,by = "Sample")
  
  #1==CRC;0==Healthy
  if(!is.null(binaryConditions)){
    HanniganGD.df = HanniganGD.df %>% 
      filter(Status%in% binaryConditions) %>% 
      mutate(Status=as.factor(if_else(Status==positiveClass,1,0)))
  }
  
  #new added5/2/19
  rownames(HanniganGD.df) = HanniganGD.df$Sample
  ####################################################3%>% 
  
  HanniganGD.df = dplyr::select(HanniganGD.df,-Sample) 
  metaData = HanniganGD.metaData %>% 
    filter(subjectID %in% rownames(HanniganGD.df))
  
  return(
    
    list(otuTable = HanniganGD.df,
         pt_MetaData = metaData,
         phyloTree = data.tree@phy_tree,
         Info = loman.eset)
  )
}