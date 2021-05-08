dispersionShift = function(tbl2,shiftedGroup,dispFact = 1,a = round(seq(0,10,by = .1),2)){
  
  sampleFromEmpDistr = function(vec,numSamples,precision = 0.0001){
    empDistr = ecdf(vec)
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
  
  
  distAtAlpha_n = function(alpha,pertb_dir,basePoint){
    ph = compositions::clo( basePoint  * compositions::clo(pertb_dir^alpha) ) 
    
    sqrt(sum( compositions::clr(ph)^2) )
  }
  
  distAtAlp = function(bp,pertDir,alp,pertPoint){
    ph = compositions::clo( bp  * compositions::clo(pertDir^alp) ) 
    pp = (compositions::clo(ph/as.numeric(bp)))
    sqrt(sum(compositions::clr(pp)^2))
  }
  
  groups = unique(tbl2[,1])
  i = which(groups==shiftedGroup)
  otherGroup = as.character(groups[-i])
  

  ### Split Data Distributions
  comp_shifted = tbl2[tbl2[,1]==shiftedGroup,-1]
  comp_fixed = tbl2[tbl2[,1]!=shiftedGroup,-1]
  
  ## centroids
  compMean.shifted = mean.acomp(acomp(comp_shifted))
  compMean.fixed = mean.acomp(acomp(comp_fixed))
  
  ### Dispersion Lengths
  ### dispersion we have for fixed distribution
  dspLen.fixed = data.frame()
  
  for(i in 1:nrow(comp_fixed)){
    ph = clo(as.numeric(comp_fixed[i,])/as.numeric(compMean.fixed))
    ph = data.frame(pointNum = i,disCentr = sqrt(sum(clr(ph)^2)))
    dspLen.fixed = rbind(dspLen.fixed,ph)
  }

  ### disoersion distibution required for shifted distribution
  newDisp = dispFact*replicate(nrow(comp_shifted),sampleFromEmpDistr(dspLen.fixed$disCentr,numSamples = 1))

  
  
  ### dispersion of the shifted distribution
  dspLen_shifted = data.frame()
  for(i in 1:nrow(comp_shifted)){
    ph = clo(as.numeric(comp_shifted[i,])/as.numeric(compMean.shifted))
    ph = data.frame(pointNum = i,disCentr = sqrt(sum(clr(ph)^2)))
    dspLen_shifted = rbind(dspLen_shifted,ph)
  }

  
  
  ########################################################
  shiftsPoints = foreach(i = 1:nrow(comp_shifted),.combine = rbind,.packages = "data.table")%dopar%{
   
    ## Define base point 
    pt = ((as.matrix(comp_shifted[i,])))
    ## compute linear shift direction
    pertb = compositions::clo(as.numeric(compMean.shifted) / as.numeric(pt))
   
    ## sample dispersion distance
    shiftDist = newDisp[i]
    
    # compute required alpha
    alpha = shiftDist/sqrt(sum(clr(pertb)^2))
    
    # ## find required alpha for sampled shift distance ****** This probably can be solved with algebra look into in the future if time permits*****
    # ## Approx solution to solving for alpha
    # simAlpha = sapply(a, function(x) distAtAlp(bp = as.numeric(compMean.shifted),pertDir = pertb,alp = x,pertPoint = pt))
    # dt = data.table(X = a,prob = simAlpha) # you'll see why val is needed in a sec
    # #setattr(dt, "sorted", "prob")  # let data.table know that w is sorted
    # setkey(dt, prob) # sorts the data
    # # binary search and "roll" to the nearest neighbour
    # # In the final expression the val column will have the you're looking for.
    # nn = shiftDist
    # xx = dt[J(nn), roll = "nearest"]
    # alpha = xx$X
    
    ## compute final shifted point at the solved alpha
    ph = compositions::clo( as.numeric(compMean.shifted)  * compositions::clo(pertb^alpha) ) 
    
    
    ph = data.frame(sampleNum = i,
                    reqAlpha = alpha,
                    actDist =  sqrt(sum(compositions::clr(compositions::clo(ph/as.numeric(compMean.shifted)))^2)),
                    requiredDist = shiftDist,
                    t(ph))
    
    ph
   
  }
  
  
  
  pointAtAlpha_n = shiftsPoints[,-4:-1]
  colnames(pointAtAlpha_n) = colnames(comp_shifted)
  
  ### New Dispersion Distr
  dspLen_after = data.frame()
  for(i in 1:nrow(pointAtAlpha_n)){
    ph = clo(as.numeric(pointAtAlpha_n[i,])/as.numeric(compMean.shifted))
    ph = data.frame(pointNum = i,disCentr = sqrt(sum(clr(ph)^2)))
    dspLen_after = rbind(dspLen_after,ph)
  }

  
  ## from real data
  baryCenter = clo(rep(1,ncol(tbl2[,-1])))
  
  tbl = rbind(data.frame(Type = otherGroup,comp_fixed),
              data.frame(Type = "dispAdjusted",pointAtAlpha_n),
              data.frame(Type = shiftedGroup,comp_shifted))
  bc = data.frame(Type ="baryCenter",t(baryCenter));colnames(bc) = colnames(tbl)
  tbl = rbind(tbl, bc)
  
  #dispersion test
  simDat = rbind(
    data.frame(Type = "dispAdjusted",pointAtAlpha_n),
    data.frame(Type = otherGroup,comp_fixed))
  
  
  
  return(list(combinedOutput = tbl,dispAdjustedData = simDat, 
         dispDistr.fixed = dspLen.fixed ,
         dispDistr.shifted = dspLen_shifted, 
         dispDistr.adj = dspLen_after   ))
}






dispersionShift_arb = function(tbl2,dispFact = 1){
  
 
  
  ## centroids
  compMean.shifted = mean.acomp(acomp(tbl2))
  

  ### dispersion of the shifted distribution
  dspLen_shifted = data.frame()
  for(i in 1:nrow(tbl2)){
    ph = clo(as.numeric(tbl2[i,])/as.numeric(compMean.shifted))
    ph = data.frame(pointNum = i,disCentr = sqrt(sum(clr(ph)^2)))
    dspLen_shifted = rbind(dspLen_shifted,ph)
  }
  
  
  #hist(dspLen.fixed$disCentr)
  d = dspLen_shifted$disCentr

  
  
  ########################################################
  shiftsPoints = foreach(i = 1:nrow(tbl2),.combine = rbind,.packages = "data.table")%dopar%{
    
    ## Define base point 
    pt = ((as.matrix(tbl2[i,])))
    ## compute linear shift direction
    pertb = compositions::clo(as.numeric(compMean.shifted) / as.numeric(pt))
    
    ## sample dispersion distance
    shiftDist = d[i]*rnorm(1,dispFact,dispFact/3)
    
    # compute required alpha
    alpha = shiftDist/sqrt(sum(clr(pertb)^2))
   
    
    ## compute final shifted point at the solved alpha
    ph = compositions::clo( as.numeric(compMean.shifted)  * compositions::clo(pertb^alpha) ) 
    
    
    ph = data.frame(sampleNum = i,
                    reqAlpha = alpha,
                    actDist =  sqrt(sum(compositions::clr(compositions::clo(ph/as.numeric(compMean.shifted)))^2)),
                    requiredDist = shiftDist,
                    t(ph))
    
    ph
    
  }
  
  
  
  pointAtAlpha_n = shiftsPoints[,-4:-1]
  colnames(pointAtAlpha_n) = colnames(tbl2)
  
  ### New Dispersion Distr
  dspLen_after = data.frame()
  for(i in 1:nrow(pointAtAlpha_n)){
    ph = clo(as.numeric(pointAtAlpha_n[i,])/as.numeric(compMean.shifted))
    ph = data.frame(pointNum = i,disCentr = sqrt(sum(clr(ph)^2)))
    dspLen_after = rbind(dspLen_after,ph)
  }
  
  

  
  #dispersion test
  simDat = rbind(
    data.frame(Type = "dispAdjusted",pointAtAlpha_n),
    data.frame(Type = "Base",tbl2))
  
  
  
  return(list(combinedOutput = simDat,
              dispAdjustedData = pointAtAlpha_n, 
             
              dispDistr.shifted = dspLen_shifted, 
              dispDistr.adj = dspLen_after   ))
}
