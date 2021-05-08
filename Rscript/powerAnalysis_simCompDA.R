rm(list=ls())
gc()

### Load Required Packages  ####
library(compositions)
library(data.table)
library(reshape2)
library(doParallel)
library(igraph)
library(caret)
library(tidyverse)
library(vegan)
library(PRROC)
library(energy)
library(Matrix)
library(ALDEx2)




clus <- parallel::makeCluster(10) 
doParallel::registerDoParallel(clus)
source(file = "C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/selectionEnergyPermutation/Functions/functions1.R")
source(file = "C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/selectionEnergyPermutation/Functions/fasterDCV_real.R")


###============================*
## Sim Data ####
###============================*
## parms
nreps  = 50
ss  =seq(25,300,by = 25)

for(sampleSize in ss){
  for(dataset in 1:4){
    powerAnalysis = data.frame()
    fname = paste0("Results/powerAnalysis/dataset",dataset,"_sampleSize",sampleSize,".csv")
    
    system.time({
      for(sd in 1:nreps){
        ###============================*
        ## exp data ####
        ###============================*
        set.seed(sd)
        switch (dataset,
                {df = read_csv("Output/16S_GEMS_microbiomeDB.csv");df = sample_n(df,size = sampleSize)},
                {df = read_csv("Output/16S_healthyBMI_microbiomeHD.csv");df = sample_n(df,size = sampleSize)},
                {df = read_csv("Output/WGS_oralCavity_HMP2012.csv");df = sample_n(df,size = sampleSize)},
                {df = read_csv("Output/WGS_controls_vataneen2016.csv");df = sample_n(df,size = sampleSize)}
        )
        
        
        
        
        ### ----------- Process data -------------
        df = data.frame(df)
        proc_dat = processCompData(tbl = df,minPrevalence = .9)
        dat = proc_dat$processedData
        y = dat[,-1]
        bool = colSums(y)==0
        y = y[,!bool]
        ## processed data.frame output
        dat = data.frame(Status = dat[,1],y)
        xx = data.frame(Status = dat[,1],fastImputeZeroes(dat[,-1]))[,-1]
        pts = simCompData(comp = xx,ptPerSample = 3,fct = 1,alFactor = 1,a = alp,modthres = .3)
        xx = sample_n(pts$simData,size = sampleSize)
        
        ## baseline shift
        pts = linearShift_arb(tbl2 = xx,a = 1,alpha_n = 1,directionType = "centroid")
        
        ### linear shift
        ### for method to detect linear shifts
        ## Linear Shift
        adjData = data.frame(Alpha = -1,Distance = 0,xx)
        a_ = seq(0,25,length.out = 10)/pts$meanDelta[2]
        new = xx
        for(i in 1:length(a_)){
          ls = linearShift_arb(tbl2 = xx,a = a_[i],alpha_n = a_[i],directionType = "centroid")
          pts = simCompData(comp =ls$linAdjustedData,ptPerSample = 3,fct = 1,alFactor = 1,a = alp,modthres = .3)
          new = data.frame(Alpha = ls$meanDelta[1],Distance = ls$meanDelta[2],sample_n(pts$simData,size = sampleSize))
          adjData = rbind(adjData,new)
        }
        
        
        
        ## compute power
        alpha = 0.05
        diffCentroid = unique(adjData$Distance)
        sampleSize = nrow(df)
        results = foreach(i = 1:length(a_),.combine = rbind)%dopar%{
          ph = adjData[adjData$Alpha %in% c(-1,a_[i]),]
          ## dispersion test
          labels = factor(ph$Alpha)
          ec  = clr(fastImputeZeroes(ph[,-2:-1]))
          d = parallelDist::parDist(as.matrix(ec))
          mod = vegan::betadisper(d,group = labels)
          vv = vegan::permutest(mod)
          
          #test equaliy of distribution
          tb = table(ph[,1])
          classes = unique(ph[,1])
          sz = c(nrow(xx) , nrow(xx) )
          et=energy::eqdist.etest(x = d,sizes = sz,distance = T,R = 1000)
          
          ## Permanova
          a.df = data.frame(Type = labels)
          pmv = vegan::adonis2(d~Type,data = a.df,permutations = 1000)
          pv = pmv$`Pr(>F)`[1]
          
          data.frame(Dist = diffCentroid[i],Alpha = a_[i],sample_Size = sampleSize,data_set = dataset,seed = sd,
                     dispPval = vv$tab$`Pr(>F)`[1],
                     p = pv, 
                     singf = if_else(pv<alpha,1,0),
                     energyTest = et$p.value,et.signif =if_else(et$p.value<alpha,1,0)
                     
          )
          
        }
        
        powerAnalysis = rbind(powerAnalysis,results)
        
        write_csv(powerAnalysis,path = fname)
        message("seed",paste0(sd,"-",fname))
      }
      
    })
  }
}


