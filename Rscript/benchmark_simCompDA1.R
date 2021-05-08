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

## selEnergy Functions/Package
source(file = "C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/selectionEnergyPermutation/Functions/functions1.R")
source(file = "C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/selectionEnergyPermutation/Functions/fasterDCV_real.R")

## simCompDa Functions.pAckage
fnames = dir(path = "Functions/")
for(f in fnames){
 source(paste0("Functions/",f) )
}


sampleSIze = 250
sd = 21

###============================*
## Sim Data ####
###============================*
## parms
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
df = read_csv("Data/16S_GEMS_microbiomeDB.csv");df = sample_n(df,size = sampleSIze)
df = read_csv("Data/16S_healthyBMI_microbiomeHD.csv");df = sample_n(df,size = sampleSIze)
df = read_csv("Data/WGS_oralCavity_HMP2012.csv");df = sample_n(df,size = sampleSIze)
df = read_csv("Data/WGS_controls_vataneen2016.csv");df = sample_n(df,size = sampleSIze)


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

visualizeData_pca(tbl = dat,sparsePercent = .9,is_relativeAbundance = T)




### ----------- Simulate Data -------------
set.seed(sd)
simulationMethod = 1
switch(simulationMethod,
       
       {  ## 1 - SimCompDA
         
         ## Stage 1
         xx = data.frame(Status = dat[,1],fastImputeZeroes(dat[,-1]))
         alp = seq(-3,3,.025)
         pts = simCompData(comp = xx[,-1],ptPerSample = 3,fct = 1,alFactor = 2,a = alp,modthres = .3);
         alp.df = pts$alphaData
         sampledPoints = pts$simData
         
         ## Stage 2
         pts = simCompData(sampledPoints,ptPerSample = 3,fct = .25,alFactor = .5,a = alp,modthres = .3);
         sampledPoints = sample_n(pts$simData,size = nrow(xx));
         
         
         ### Combined Empirical and Simulated Data Table
         tbl = rbind(data.frame(Type = "Empirical",(dat[,-1])),
                     data.frame(Type = "Simulated",(sampledPoints)))
         
        
         
         labels = tbl[,1]
         
       },
       {  ## 2 - dirMult
         
         dm_fit = DirichletMultinomial::dmn(as.matrix(dat[,-1]),k = 1)
         est = dm_fit@fit$Estimate[,1]
         
         ## fit model
         # f = MGLM::MGLMfit(data = dat[,-1],dist = "DM")
         # est = f@estimate
         
         ## simulate from fit
         numSamps = 1 * nrow(dat)
         totalCounts = rowSums(y)
         Ysim <- MGLM::rdirmn(n = numSamps, 
                              size = sample(totalCounts,size = numSamps,replace = T),
                              alpha = est)
         colnames(Ysim) = colnames(y)
         tbl = rbind(data.frame(Type = "Empirical",clo(y)),
                     data.frame(Type="Simulated",clo(Ysim)))
         #Test statistc
         labels = tbl[,1]
         
       },
       {  ## 3 - zinb
         
         ft.zinb =zinbwave::zinbFit(t(dat[,-1]),zeroinflation = T)
         nreps = 1
         Ysim1  = data.frame()
         for(i in 1:nreps){
           Ysim = t(zinbwave::zinbSim(ft.zinb)$counts)
           colnames(Ysim) = colnames(y)
           Ysim1 = rbind(Ysim1,Ysim)
         }  
         colnames(Ysim) = colnames(dat[,-1])
         tbl = rbind(data.frame(Type = "Empirical",clo(dat[,-1])),data.frame(Type="Simulated",clo(Ysim1)))
         labels = simDat[,1]
         
       },
       
       { ## 4 - ALDEX2
         Ysim = data.frame()
         nreps = 1
         for(j in 1:nreps){
           x <- aldex.clr(t(dat[,-1]) , dat[,1], mc.samples=1, denom="all", verbose=F,)
           ald_dat = x@analysisData
           ph.sim = data.frame()
           for(i in 1:length(ald_dat)){
             ph = t(ald_dat[[i]])
             ph.sim = rbind(ph.sim,ph)
           }
           ph.sim = data.frame(Type = "Simulated",compositions::clrInv(ph.sim))
           Ysim = rbind(Ysim,ph.sim)
         }
         tbl = rbind(data.frame(Type = "Empirical",fastImputeZeroes(y) ),Ysim)
         labels = tbl[,1]
       }
       
)





### ----------- Visulalize Data -------------

## Unweighted PCA
pc = prcomp(clr(fastImputeZeroes(tbl[,-1])),center = F,scale. = T)
coords.df = data.frame(Type = tbl[,1],pc$x)
ggplot(coords.df,aes(PC1,PC2,))+
  geom_point(aes(shape = Type,col= Type,fill = Type,alpha = Type),size = 3)+
  scale_shape_manual(values = c(16,22))+
  scale_alpha_manual(values = c(1,.5))+
  #geom_text()+
  #scale_color_viridis_c(option = "E")+
  theme_bw()+
  #stat_density_2d(n = 200)+
  ggtitle("Zero Inflated Negative Binomial Simulation",subtitle = "WGS Microbiome Data")

## Weighted CLR PCA
ec  = easyCODA::CLR(fastImputeZeroes(tbl[,-1]),weight = T)
pc = easyCODA::PCA(ec)
pc_ = pc$rowcoord
coords.df = data.frame(Type = tbl[,1],pc_)
#orig_comp = coords%*%solve(X.evct)
ggplot(coords.df,aes(X1,X2,col= Type,Fill = Type))+
  geom_point(aes(shape = Type,col= Type,fill = Type,alpha = Type),size = 3)+
  scale_shape_manual(values = c(16,22))+
  scale_alpha_manual(values = c(1,.5))+
  #geom_text()+
  #scale_color_viridis_c(option = "E")+
  theme_bw()+
  #stat_density_2d(n = 200)+
  ggtitle("Zero Inflated Negative Binomial Simulation",subtitle = "WGS Microbiome Data")






### ----------- Compare Distribution Chararceristics  -------------

## compute CLR transform
ph = data.frame(Type = tbl[,1],ilr(fastImputeZeroes(tbl[,-1])))

## Correlatation between mean vectors
truCov = colMeans(ph[ph[,1] == "Empirical",-1])
simCov = colMeans(ph[ph[,1] != "Empirical",-1])
lepage.test(truCov,simCov)
cor(truCov,simCov,method = "spearman")
plot(truCov,simCov)

# ## difference in correlation
truCov = cor(ph[ph$Type=="Empirical",-1],method = "spearman")
simCov = cor(ph[ph$Type!="Empirical",-1],method = "spearman")
cor(as.vector(truCov[lower.tri(truCov)]) , as.vector(simCov[lower.tri(simCov)]) ,method = "spearman")
plot(as.vector(truCov[lower.tri(truCov)]) , as.vector(simCov[lower.tri(simCov)]) )


## simultan diff in cov or mean
{labels = tbl[,1]
  classes = unique(tbl[,1])
  alr.simdat = compositions::alr(fastImputeZeroes(tbl[,-1]))
  m1 =  colMeans(alr.simdat[labels==classes[1],])
  m2 = colMeans(alr.simdat[labels==classes[2],])
  s1 = cov(alr.simdat[labels==classes[1],])
  s2 = cov(alr.simdat[labels==classes[2],])
  N1 = nrow(alr.simdat[labels==classes[1],])
  N2 = nrow(alr.simdat[labels==classes[2],])
  sp =  (N1*s1 +N2*s2) / (N1+N2)
  mc = (N1*m1 +N2*m2) / (N1+N2)
  sc = sp + ((N1+N2)^-2) * N1*N2*(m1-m2)%*%t(m1-m2)
  tsat = ((N1*log(norm(sc,type = "F")/norm(s1,type = "F"))) + (N2 * log(norm(sc,type = "F")/norm(s2,type = "F"))))^2
  ((N1*log(norm(sc,type = "2")/norm(s1,type = "2"))) + (N2 * log(norm(sc,type = "2")/norm(s2,type = "2"))))
  # difference in covariance 
  l2_tstat.cov =  ((N1*log(norm(sp,type = "F")/norm(s1,type = "F"))) + (N2 * log(norm(sp,type = "F")/norm(s2,type = "F"))))^2
  nreps = 100
  nulldistr_tstat = foreach(x = 1:nreps,.combine = rbind)%dopar%{
    labels = sample(labels)
    alr.simdat = compositions::alr(fastImputeZeroes(tbl[,-1]))
    m1 =  colMeans(alr.simdat[labels==classes[1],])
    m2 = colMeans(alr.simdat[labels==classes[2],])
    s1 = cov(alr.simdat[labels==classes[1],])
    s2 = cov(alr.simdat[labels==classes[2],])
    N1 = nrow(alr.simdat[labels==classes[1],])
    N2 = nrow(alr.simdat[labels==classes[2],])
    sp =  (N1*s1 +N2*s2) / (N1+N2)
    mc = (N1*m1 +N2*m2) / (N1+N2)
    sc = sp + ((N1+N2)^-2) * N1*N2*(m1-m2)%*%t(m1-m2)
    ((N1*log(norm(sc,type = "F")/norm(s1,type = "F"))) + (N2 * log(norm(sc,type = "F")/norm(s2,type = "F"))))^2
  }
  critical_value.u = quantile(nulldistr_tstat,probs = 0.975)
  critical_value.l = quantile(nulldistr_tstat,probs = 0.025)
  #reject null if tstat > critical value
  
  (sum(tsat<nulldistr_tstat)+1) / (nreps+1)
  hist(nulldistr_tstat,breaks = 50)
  abline(v = tsat,col = "red",lty = "dashed")
}## Lattice level 2
(sum(tsat<nulldistr_tstat)+1) / (nreps+1)

## Covariance diff
{nulldistr_l2_tstat = foreach(x = 1:nreps,.combine = rbind)%dopar%{
  labels = sample(labels)
  alr.simdat =compositions::alr(fastImputeZeroes(tbl[,-1]))
  m1 =  colMeans(alr.simdat[labels==classes[1],])
  m2 = colMeans(alr.simdat[labels==classes[2],])
  s1 = cov(alr.simdat[labels==classes[1],])
  s2 = cov(alr.simdat[labels==classes[2],])
  N1 = nrow(alr.simdat[labels==classes[1],])
  N2 = nrow(alr.simdat[labels==classes[2],])
  sp =  (N1*s1 +N2*s2) / (N1+N2)
  mc = (N1*m1 +N2*m2) / (N1+N2)
  sc = sp + ((N1+N2)^-2) * N1*N2*(m1-m2)%*%t(m1-m2)
  ((N1*log(norm(sp,type = "F")/norm(s1,type = "F"))) + (N2 * log(norm(sp,type = "F")/norm(s2,type = "F"))))^2
}
  critical_value = quantile(nulldistr_l2_tstat,probs = 0.95)
  critical_value.u = quantile(nulldistr_l2_tstat,probs = 0.975)
  critical_value.l = quantile(nulldistr_l2_tstat,probs = 0.025)
  
  #reject null if tstat > critical value; i.e the covaraince between the two distribution is not equal ; Ha is that they are equal
  #Fail to reject null means distr are different
  # Reject null means distr are the same
  (sum(l2_tstat.cov<nulldistr_l2_tstat)+1) / (nreps+1)
  hist(nulldistr_l2_tstat,breaks = 30)
  abline(v = l2_tstat.cov,col = "red",lty = "dashed")
}
(sum(l2_tstat.cov<nulldistr_l2_tstat)+1) / (nreps+1)
# # difference in mean using two-sample location scale problem.
lp.test = lepage.test(m1,m2)




### ----------- Beta Diveristy -------------

## dispersion test
labels = tbl[,1]
ec  = easyCODA::CLR(fastImputeZeroes(tbl[,-1]),weight = F)
d = parallelDist::parDist(as.matrix(ec$LR))
mod = vegan::betadisper(d,group = labels)
anova(mod)
permutest(mod)
mod.HSD <- TukeyHSD(mod)
plot(mod,)
boxplot(mod)
plot(TukeyHSD(mod))
plot(mod, ellipse = F, hull = FALSE,label = F) # 1 sd data ellipse

## Permanova
a.df = data.frame(Type = labels)
pmv = adonis2(d~Type,data = a.df,permutations = 100)
(pmv)




### ----------- Compare Alpha Diversity Characteristics -------------

## Compte Alpha DIveristy
classes = unique(tbl[,1])
data_ = fastImputeZeroes(tbl[,-1])
shannon <- vegan::diversity(data_, "shannon")
simp =  vegan::diversity(data_, "simpson")
invsimp <-  vegan::diversity(data_, "inv")
diversity.df = data.frame(Group = tbl[,1],
                          H = shannon,
                          Simpson = simp,
                          invSimpson = invsimp
)
diversityIndex = "Simpson"
div.gathered = diversity.df %>% 
  gather(key = "Diversity",value = "Index",2:ncol(diversity.df)) %>% 
  filter(Diversity == diversityIndex)
my_comparisons <- list( c("Simulated","Empirical"))

## Compare using wilcox test 
wt = wilcox.test(x = div.gathered[div.gathered[,1]==classes[1],3],
            y =  div.gathered[div.gathered[,1]==classes[2],3]
)
wt$p.value

ks.test(x = div.gathered[div.gathered[,1]==classes[1],3],
        y =  div.gathered[div.gathered[,1]==classes[2],3])

## Plot diveristy
ggplot(div.gathered,aes(Group,Index,col = Group))+
  geom_boxplot(outlier.shape = NA)+
  scale_color_manual(values = c("Red","Blue"))+
  geom_jitter(aes(fill = Group),width = .1,size = 3,alpha = .1,
              #shape = 21,col = "black"
  )+
  scale_fill_manual(values = c("Red","Blue"))+
  theme_classic()+
  theme(legend.position = "none")+
  ylab("Shannon Entropy")


  
### ----------- Machine Learning Comparision -------------
## ML model
#xx = data.frame(Status = tbl[,1],(fastImputeZeroes(tbl[,-1])))
xx = data.frame(Status = tbl[,1],clr(fastImputeZeroes(tbl[,-1])))
#xx = data.frame(Status =tbl[,1],ec$LR)

trainData = xx
testData = xx
yt = (trainData[,1])
md = trainML_Models(trainLRs = trainData[,-1],ytrain = yt,testLRs = testData[,-1],
                    y_test = testData[,1],testIDs = data.frame(ID = 1:nrow(testData)),
                    models = c("lda","knn","gbm","svmLinear","nnet"),
                    numRepeats = 1,
                    numFolds = 10)

md$performance
mean(md$performance$TrainAUC)

