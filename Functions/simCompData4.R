
simCompData = function(comp,labels = NULL, coarseGrained = T,a = seq(0,10,by = .05),k = round(sqrt(nrow(comp))),  nreps = 100,modthres = .4,ndirections = 1000,fct = 1.5, alFactor = 1.5,ptPerSample = 3){
  
  distAtAlp = function(bp,pertDir,alp,pertPoint){
    ph = clo( bp  * clo(pertDir^alp) ) 
    pp = (clo(ph/as.numeric(bp)))
    sqrt(sum(clr(pp)^2))
  }
  

    ## compute pairwise distance 
    dd = parallelDist::parallelDist(clr(comp))
    ## Compyte KNN Topology and communiyt structure
    g = knn_graph(as.matrix(dd),K = k,sim_ = T)
    #Extract KNN Adj Matrix
    adjacency_matrix <- igraph::as_adjacency_matrix(graph = g$Graph,sparse = F,
                                                    attr = "weight"
    )
    ## Define Community Structure
    partition <- leiden::leiden(object = adjacency_matrix)
    cls = data.frame(membership = partition)
    testMod = igraph::modularity(g$Graph,partition)
    plot(g$Graph,vertex.color =factor(cls$membership),vertex.label = NA,vertex.size = 5,edge.width = .5)

    
    ## if significant use community structure
    if(testMod<modthres){
      message("Not signif.", "; modularity = ",round(testMod,3))
      cls = data.frame(membership = rep(1,nrow(comp)))
      comIDs = 1
    }else{
      message("Signif.", "; modularity = ",round(testMod,3))
      ss = table(cls$membership)
      
      ## combine groups with single member
      ii = which(as.numeric(ss)==1)
      if(!is_empty(ii)){
        if(length(ii)==1){
          comb = as.numeric(names(ss)[ii])
          cls$membership[cls$membership %in%comb] = max(unique(cls$membership))
        }else{
          comb = as.numeric(names(ss)[ii])
          cls$membership[cls$membership %in%comb] = max(unique(cls$membership))+1  
        }
      }
    }
      
    comIDs = unique(cls$membership)
    comIDs = sort(comIDs)
    
  ## Data Frames
  shiftsPoints = data.frame()
  reqAlpha = data.frame()
  allDat = data.frame()
  
  message("Compute Augmentation....")
  for(cid in comIDs){
    
    ## select community points
    phd = data.frame(comp[cls$membership==cid,])
    phd.labels = labels[cls$membership==cid]
    est = compositions::fitDirichlet(acomp(phd))
    alpha_ = est$alpha
    simDir_vec = compositions::rDirichlet.acomp(ndirections,alpha = alFactor*alpha_)#c(mean_comp)*200)
    
    compMean.fixed = mean.acomp(acomp(phd))
    dspLen.fixed = data.frame()
    for(i in 1:nrow(phd)){
      ph = clo(as.numeric(phd[i,])/as.numeric(compMean.fixed))
      ph = data.frame(pointNum = i,disCentr = sqrt(sum(clr(ph)^2)))
      dspLen.fixed = rbind(dspLen.fixed,ph)
    }
    d = dspLen.fixed$disCentr
    cv = sd(d) / mean(d)
    mn_d = mean(d*fct*cv); sd_d = sd(d*fct*cv)

    for(r in 1:ptPerSample){
      
      ad = foreach(i = 1:nrow(phd),.combine = rbind,.packages = c("tidyverse","compositions","data.table"))%dopar%{
        
        ## seelct point
        pt = ((as.matrix(phd[i,])))
        ## randomly sample a direction
        rx = sample(x = 1:nrow(simDir_vec),size = 1)
        randDir = simDir_vec[rx,]
        ## define pertubation from point to new direction
        pertb = clo(as.numeric(pt) / as.numeric(randDir)) ## from point -> rand dir
        ## sample shift distance
        shiftDist = rnorm(1,mn_d,sd_d)
        
        ## compute the required alpha
        alpha = shiftDist/sqrt(sum(clr(pertb)^2))
        
        ## Get alpha
        if(runif(1)>.5){
          alpha = alpha*-1
        }else{
          alpha = alpha
        }
        
        
        ## define new point in random direction of length shift dist @ alpha n
        ph = clo( pt  * clo(pertb^alpha) ) 
        
        
        if(is.null(labels)){
          alphaData_debug = data.frame(rep = r,sampleNum = i,
                                       comm = cid,
                                       actAlpha = alpha,
                                       act_shiftDistance = sqrt(sum(clr(clo(ph/pt))^2)),
                                       diff = abs(sqrt(sum(clr(clo(ph/pt))^2)) - shiftDist),
                                       reqShift = shiftDist)
        }else{
          alphaData_debug = data.frame(classLabel = phd.labels[i],rep = r,sampleNum = i,
                                       comm = cid,
                                       actAlpha = alpha,
                                       act_shiftDistance = sqrt(sum(clr(clo(ph/pt))^2)),
                                       diff = abs(sqrt(sum(clr(clo(ph/pt))^2)) - shiftDist),
                                       reqShift = shiftDist)
        }
       
        
      
        
        ## Combine Data 
        cbind(alphaData_debug,ph)

      }
      
      allDat = rbind(allDat,ad)
      message("Rep = ",r ," of ",ptPerSample,"; Comm =  ",cid," of ",length(comIDs))
    }
    
  }
  

  
 if(is.null(labels)){
   classLabel = data.frame()
   reqAlpha = allDat[,1:7]
   shiftsPoints = allDat[,-7:-1]
 }else{
   classLabel = allDat[,1]
   allDat = allDat[,-1]
   reqAlpha = allDat[,1:7]
   shiftsPoints = allDat[,-7:-1]
 }

  return(list(simData = shiftsPoints, alphaData = reqAlpha,ClassLabels = classLabel ))
}

