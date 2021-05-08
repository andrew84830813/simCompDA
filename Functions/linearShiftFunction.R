vecNorm = function(namedVec){
  c = combinat::combn2(names(namedVec))
  el  = data.frame(c,ff = paste0(c[,1],"___",c[,2]))
  lrs = getLogRatios(data.frame(t(namedVec)),ratioList = el)
  sqrt(sum(lrs^2)/length(lrs))
}

vecDist = function(v1,v2){
  c = combinat::combn2(names(v1))
  el  = data.frame(c,ff = paste0(c[,1],"___",c[,2]))
  
  v1 = getLogRatios(data.frame((v1)),ratioList = el)
  
  v2 = getLogRatios(data.frame((v2)),ratioList = el)
  
  sqrt(sum((v1-v2)^2)/length(v1))
}

linearShift = function(tbl2,shiftedGroup,pertubationDirection = NULL,a =1,alpha_n = 1){
  
  sampleFromEmpDistr = function(vec,numSamples,precision = 0.0001){
    empDistr = ecdf(dspLen$disCentr)
    x = seq(min(vec),max(vec),precision)
    ## sample from empirical distribution
    dt = data.table(X = x,prob = empDistr(x)) # you'll see why val is needed in a sec
    setattr(dt, "sorted", "prob")  # let data.table know that w is sorted
    setkey(dt, prob) # sorts the data
    # binary search and "roll" to the nearest neighbour
    # In the final expression the val column will have the you're looking for.
    nn = runif(numSamples)
    xx = dt[J(nn), roll = "nearest"]
    xx$X
  }
  
  
  fn.dist_wAlpha=function(base,direction,alpha = 1){
    comp = compositions::clo(base/compositions::clo(direction^alpha))
    z = as.matrix(comp)
    s=lapply(1:ncol(z), function(x) log(z/z[,x]))
    s=sapply(1:ncol(z),function(x) s[[x]]^2)
    sqrt(sum(s[lower.tri(s)]))
  }
  
  groups = unique(tbl2[,1])
  i = which(groups==shiftedGroup)
  otherGroup = as.character(groups[-i])
  message(otherGroup)
  
  distAtAlpha_n = function(alpha,pertb_dir,basePoint){
    ph = clo( basePoint  * clo(pertb_dir^alpha) ) 
    
    sqrt(sum(clr(ph)^2))
  }
  
  distAtAlp = function(bp,pertDir,alp,pertPoint){
    ph = clo( bp  * clo(pertDir^alp) ) 
    pp = (clo(ph/as.numeric(bp)))
    sqrt(sum(clr(pp)^2))
  }
  
  
  ### Split Data Distributions
  comp_shifted = tbl2[tbl2[,1]==shiftedGroup,-1]
  comp_fixed = tbl2[tbl2[,1]!=shiftedGroup,-1]
  
  ## centroids
  compMean.shifted = mean.acomp(acomp(comp_shifted))
  compMean.fixed = mean.acomp(acomp(comp_fixed))
  
  ## Direction to mean 
  if(is.null(pertubationDirection)){
    pertb = as.vector(clo(data.frame(t(compMean.fixed))/data.frame(t(compMean.shifted))))
  }
  
  ########################################################
  #Alpha Range

  shiftsPoints = data.frame()
  
  shiftsPoints = foreach(i = 1:nrow(comp_shifted),.combine = rbind)%dopar%{
    c.df_pertb = data.frame()
    for(iA in a){
      ph = compositions::clo( ((as.matrix(comp_shifted[i,])))  * compositions::clo(pertb^iA) ) 
      c.df_pertb = rbind(c.df_pertb,ph)
    }
    data.frame(sampleNum = i,Alpha = a,c.df_pertb)
  }
  
  ## select point cloud relative to specific alpha
  pointAtAlpha_n = shiftsPoints[shiftsPoints$Alpha==alpha_n,-2:-1]

  ## from real data
  baryCenter = clo(rep(1,ncol(dat[,-1])))
  
  ## combine all output
  tbl = rbind(data.frame(Type = shiftedGroup,comp_shifted),
              data.frame(Type = "linearAdjusted",pointAtAlpha_n),
              data.frame(Type = otherGroup,comp_fixed))
  # bc = data.frame(Type ="baryCenter",t(baryCenter));colnames(bc) = colnames(tbl)
  # tbl = rbind(tbl, bc)
  # 
  
  simDat = rbind(
    data.frame(Type = "linearAdjusted",pointAtAlpha_n),
    data.frame(Type = otherGroup,comp_fixed)) 
  
  return(list(combinedOutput = tbl,linAdjustedData = simDat, allShiftData = shiftsPoints) 
         )
}




linearShift_arb = function(tbl2,pertubationDirection = NULL,directionType = "centroid",evComp = 1,a = 1,alpha_n = 1){
  
  ## centroids
  compMean.shifted = mean.acomp(acomp(tbl2))
  
  ## Direction to mean 
  if(is.null(pertubationDirection)){
    
    if(directionType == "centroid" || directionType== "Centroid" || directionType=="C"){
      pertb = compMean.shifted  
    }else if(directionType =="EV" || directionType =="eigen"){
      c1 = eigenDecomp.CLR(tbl2,weighted = F)
      ev.c1  = c1$eigenVector
      c1.evDir = clrInv(ev.c1[,evComp])
      ## compute pertubation
      pertb = c1.evDir
    }
    
  }
  
  
 
  ########################################################
  #Alpha Range
  
  shiftsPoints = foreach(i = 1:nrow(tbl2),.combine = rbind)%dopar%{
    c.df_pertb = data.frame()
    for(iA in a){
      ph = compositions::clo( ((as.matrix(tbl2[i,])))  * compositions::clo(pertb^iA) ) 
      c.df_pertb = rbind(c.df_pertb,ph)
    }
    data.frame(sampleNum = i,Alpha = a,c.df_pertb)
  }
  
  ## select point cloud relative to specific alpha
  pointAtAlpha_n = shiftsPoints[shiftsPoints$Alpha==alpha_n,-2:-1]
  

  
  simDat = rbind(
    data.frame(Type = "linearAdjusted",pointAtAlpha_n),
    data.frame(Type = "Base",tbl2)) 
  
  
  ## centroids
  compMean.new = mean.acomp(acomp(pointAtAlpha_n))
  
  ## distance between centroids
  
  ds = vecDist(data.frame(t(compMean.new)),data.frame(t(compMean.shifted))) #sqrt(sum( ( as.numeric(clr(compMean.new)) - as.numeric(clr(compMean.shifted)) )^2   ))
  
  return(list(combinedOutput = simDat,
              linAdjustedData = pointAtAlpha_n, 
              allShiftData = shiftsPoints,shiftVec = ( as.numeric(clr(compMean.new)) - as.numeric(clr(compMean.shifted)) )^2 ,
              meanDelta = c(alpha = a,distance=ds)) 
  )
}




