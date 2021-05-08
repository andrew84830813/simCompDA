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
sampleSIze = 250
sd = 210

mnCounts = 1e7 # WGS
mnCounts = 1e3 # 16S
dims = 100
sizeParm = sample(1:4,size = 1)
libSize = rnbinom(n = sampleSIze,size = sizeParm,mu = mnCounts);hist(libSize)

## sim distr
df = simDirDistr(n1 = sampleSIze,dms_ = dims,seed = sd,scale = 1/log(.5*dims))
df =  simAddLogNormal(n1 = sampleSIze,dms_ = dims,seed = sd,varParm = log(dims))
df = sampleAddLogisticT(n = sampleSIze,dms = dims,seed = sd,df_ = dims)

# convert to counts
xx = round(sweep(df[,-1],MARGIN = 1,STATS = libSize,FUN = "*"))
df = data.frame(Status = df[,1],xx)

###============================*
## exp data ####
###============================*
df = read_csv("Output/16S_GEMS_microbiomeDB.csv");df = sample_n(df,size = sampleSIze)
df = read_csv("Output/16S_healthyBMI_microbiomeHD.csv");df = sample_n(df,size = sampleSIze)
df = read_csv("Output/WGS_oralCavity_HMP2012.csv");df = sample_n(df,size = sampleSIze)
df = read_csv("Output/WGS_controls_vataneen2016.csv");df = sample_n(df,size = sampleSIze)


hist(rowSums(df[,-1]))


### ----------- Process data -------------

df = data.frame(df)
proc_dat = processCompData(tbl = df,minPrevalence = .9)
dat = proc_dat$processedData
y = dat[,-1]
bool = colSums(y)==0
y = y[,!bool]
## processed data.frame output
dat = data.frame(Status = dat[,1],y)

#visualizeData_pca(tbl = dat,sparsePercent = .9,is_relativeAbundance = T)
xx = data.frame(Status = dat[,1],fastImputeZeroes(dat[,-1]))[,-1]

### Summary
#comMedian = as.numeric(Gmedian::Gmedian(clr(xx[,-1])))
compMean = mean.acomp(acomp(xx))
compMean.clr  = as.numeric(clr(compMean))
dotProduct = c()
for(i in 1:nrow(xx)){
  ph = as.numeric(clr(xx[i,]))
  dotProduct[i] = sum(ph*compMean.clr)
}
hist(dotProduct)
nm = names(compMean)



## dispersion
# ph = xx
# ph = rbind(data.frame(t(compMean)),ph)
# ph = calcLogRatio(data.frame(Status = "pos",ph))[,-1]
# compmean = as.numeric(ph[1,])
# ph = ph[-1,]
# dspLen_after = data.frame()
# for(i in 1:nrow(ph)){
#   x = sqrt( sum((as.numeric(ph[i,]) - as.numeric(compmean))^2) / length(ph) )
#   x = data.frame(pointNum = i,disCentr = x)
#   dspLen_after = rbind(dspLen_after,x)
# }
# hist(dspLen_after$disCentr)


dspLen_after1 = data.frame()
for(i in 1:nrow(xx)){
  ph = clo(as.numeric(xx[i,-1])/as.numeric(compMean))
  ph = data.frame(pointNum = i,disCentr = sqrt(sum(clr(ph)^2)))
  dspLen_after1 = rbind(dspLen_after1,ph)
}
hist(dspLen_after1$disCentr)


# dspLen_after1 = data.frame()
# for(i in 1:nrow(xx)){
#   x = clo(as.numeric(xx[i,])/as.numeric(compMean))
#   names(x)=nm
#   x = data.frame(pointNum = i,disCentr = vecNorm(x))
#   dspLen_after1 = rbind(dspLen_after1,x)
# }
# hist(dspLen_after1$disCentr)

# 
# dspLen_after1 = data.frame()
# for(i in 1:nrow(xx)){
#   v1 = xx[i,]
#   names(v1) = nm  
#   
#   v2 = data.frame(t(compMean))
#   names(v2) = nm
#   
#   x = data.frame(pointNum = i,disCentr = vecDist(v1,v2))
#   dspLen_after1 = rbind(dspLen_after1,x)
# }
# hist(dspLen_after1$disCentr)


#centroid relative location
vecNorm(compMean)
x = compMean

v = rep(1e-11,length(compMean))
v[1] = 1
names(v)  = nm
vecNorm(v)
vecNorm(compMean) / vecNorm(v)




### for method to detect linear shifts
numOut = 4

## Linear Shift
adjData = data.frame(Alpha = 0,Distance = 0,xx)
new = xx
a_ = seq(0,2,length.out = numOut)

for(i in 2:length(a_)){
  ls = linearShift_arb(tbl2 = xx,a = a_[i],alpha_n = a_[i],directionType = "centroid")
  pts = simCompData(comp =ls$linAdjustedData,ptPerSample = 3,fct = 1,alFactor = 1,a = alp,modthres = .3)
  new = data.frame(Alpha = ls$meanDelta[1],Distance = ls$meanDelta[2],sample_n(pts$simData,size = nrow(xx)))
  adjData = rbind(adjData,new)
}
# plot(acomp(adjData[,-2:-1]),col = factor(adjData$Alpha))

## laplcian eigen map
library(dimRed)
xx1 <- dimRedData(data.frame(ilr(adjData[,-2:-1])))
##isomap
leim <- dimRed::Isomap()
parms = leim@stdpars
parms$knn = 25#round(sqrt(nrow(xx)))
parms$get_geod = T
emb <- leim@fun(xx1, parms)
# ### visual
coords.df = data.frame(Alpha = adjData[,1],Distance = adjData[,2],X1 = emb@data@data[,1],X2 = emb@data@data[,2])
coords.df$Alpha = factor(coords.df$Alpha)
ggplot(coords.df,aes(X1,X2,col = Distance))+
  geom_point(alpha = .95,size = 4)+
  scale_shape_manual(values = c(15,16,17,18))+
  #geom_text()+
  #scale_color_brewer(palette = "Set2")+
  #scale_color_viridis_d(option = "D",end = .5)+
  scale_color_viridis_c(option = "D",end = .75,direction = -1)+
  # scale_fill_distiller(palette = "PuBuGn")+
  # scale_color_distiller(palette = "PuBuGn")+
  theme_bw()+
  #stat_ellipse()+
  theme(legend.position = "top")



### dispersion shift
adjData = data.frame(dispFact = 1,xx[,-1])
new = xx
dispF = seq(1,.25,length.out = numOut)

for(i in 2:length(dispF)){
  ## dispersion
  ds = dispersionShift_arb(tbl2 = xx[,-1],dispFact = dispF[i])
  new = data.frame(dispFact = dispF[i],ds$dispAdjustedData)
  adjData = rbind(adjData,new)
}

#Visualize
# pc = prcomp((adjData[,-1]),center = T)
# coords.df = data.frame(dispFact = adjData$dispFact,pc$x)
# ggplot(coords.df,aes(PC1,PC2,col = dispFact,shape = factor(dispFact)))+
#   geom_point(alpha = .95,size = 4)+
#   scale_shape_manual(values = c(15,16,17,18))+
#   #geom_text()+
#   #scale_color_brewer(palette = "Set2")+
#   #scale_color_viridis_d(option = "D",end = .5)+
#   scale_color_viridis_c(option = "D",end = .75,direction = -1)+
#   # scale_fill_distiller(palette = "PuBuGn")+
#   # scale_color_distiller(palette = "PuBuGn")+
#   theme_bw()+
#   #stat_ellipse()+
#   theme(legend.position = "top",panel.grid = element_blank())

## dispersion
xx1 <- dimRedData(data.frame(ilr(adjData[,-1])))
leim <- dimRed::Isomap()
parms = leim@stdpars
parms$knn = 15#round(sqrt(nrow(xx)))
parms$get_geod = T
emb <- leim@fun(xx1, parms)
# ### visual
coords.df = data.frame(dispFact = adjData$dispFact,X1 = emb@data@data[,1],X2 = emb@data@data[,2])
#coords.df$Alpha = factor(coords.df$Alpha)
ggplot(coords.df,aes(X1,X2,col = dispFact))+
  geom_point(alpha = .95,size = 4)+
  scale_shape_manual(values = c(15,16,17,18))+
  #geom_text()+
  #scale_color_brewer(palette = "Set2")+
  #scale_color_viridis_d(option = "D",end = .5)+
  scale_color_viridis_c(option = "D",end = .75,direction = -1)+
  # scale_fill_distiller(palette = "PuBuGn")+
  # scale_color_distiller(palette = "PuBuGn")+
  theme_bw()+
  #stat_ellipse()+
  theme(legend.position = "top",panel.grid = element_blank())







### dispersion shift
adjData = data.frame(dispFact = 1,Alpha = 0,xx[,-1])

for(i in 2:length(dispF)){
  ls = linearShift_arb(tbl2 = xx[,-1],a = a_[i],alpha_n =  a_[i],directionType = "centroid")
  ## linear shft
  pts = simCompData(comp =ls$linAdjustedData,ptPerSample = 3,fct = 1,alFactor = 1,a = alp,modthres = .3)
  new = data.frame(dispFact = dispF[i],Alpha = a_[i],sample_n(pts$simData,size = nrow(xx)))
  ## dispersion
  ds = dispersionShift_arb(tbl2 = new[,-2:-1],dispFact = dispF[i])
  new = data.frame(dispFact = dispF[i],Alpha = a_[i],sample_n(ds$dispAdjustedData,size = nrow(xx)))
  adjData = rbind(adjData,new)
}


 xx1 <- dimRedData(data.frame(ilr(adjData[,-2:-1])))
# ##isomap
leim <- dimRed::Isomap()
parms = leim@stdpars
parms$knn = 25#round(sqrt(nrow(xx)))
parms$get_geod = T
emb <- leim@fun(xx1, parms)
### laplacian
# leim <- LaplacianEigenmaps()
# parms = leim@stdpars
# parms$knn = 15#round(sqrt(nrow(xx)))
# emb <- leim@fun(xx1, parms)
# # ### visual
coords.df = data.frame(dispFact = adjData$dispFact,Alpha = round(adjData$Alpha,3),X1 = emb@data@data[,1],X2 = emb@data@data[,2])
ggplot(coords.df,aes(X1,X2,col = dispFact,shape = factor(Alpha)))+
  geom_point(alpha = .95,size = 4)+
  scale_shape_manual(values = c(15,16,17,18))+
  #geom_text()+
  #scale_color_brewer(palette = "Set2")+
  #scale_color_viridis_d(option = "D",end = .5)+
  scale_color_viridis_c(option = "D",end = .75,direction = -1)+
  # scale_fill_distiller(palette = "PuBuGn")+
  # scale_color_distiller(palette = "PuBuGn")+
  theme_bw()+
  #stat_ellipse()+
  theme(legend.position = "top",panel.grid = element_blank())
# 
# 





## Dispersion Shift
ds = dispersionShift_arb(tbl2 = xx[,-1],dispFact = 1.5)
tbl = ds$combinedOutput
