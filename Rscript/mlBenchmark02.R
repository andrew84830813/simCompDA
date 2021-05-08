Rcpp::sourceCpp("C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/selectionEnergyPermutation/Functions/cpp_ksTest3.cpp")
source(file = "C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/selectionEnergyPermutation/Functions/functions1.R")
source(file = "C:/Users/andrew84/Documents/MicrobiomeProject/Data/R_Projects/Test/selectionEnergyPermutation/Functions/fasterDCV_real.R")

mean_ = function(x){
  mean(x,na.rm = T)
}

## Data Simulatinos
df = read_csv("Output/16S_IBS-unnbalanced_microbiomeHD.csv")
df = read_csv("Output/WGS_CDI_VINCENT2016.csv")
df= read_csv("Output/16S_IBDmulticlass_ihmp.csv")
df = read_csv("Output/")


### Process Data ####
## input Data
dat = data.frame(df)
## Process Sparisty
proc_dat = processCompData(tbl = dat,minPrevalence = .9)
dat = proc_dat$processedData
y = dat[,-1]
bool = colSums(y)==0
y = y[,!bool]
## processed data.frame output
dat = data.frame(Status = dat[,1],y)
dat = data.frame(Status = dat[,1],fastImputeZeroes(dat[,-1]))
minClass = proc_dat$minClss
majClass = proc_dat$majClass


## dat aug parms
balanceClasses = F; pts.sam = 4; numSamples = 2*nrow(dat)
perRows = .4
modularityThreshold = .25
percentMaj = .25




## Splilt Data Into K-fold
Outcome = colnames(dat)[1]
augData = T
nfolds = 5
sd = 109
allData = kfoldDataPartition(df = dat,
                             kfold = nfolds,
                             permuteLabel = F,
                             seed = sd)

sensOptPerformance = data.frame()

for(i in 1:nfolds){
  set.seed(i)
  
  td = allData[[i]]$xtrain_combinedFolds
  td.id = allData[[i]]$xtrain_IDs
  ts = allData[[i]]$xtest_kthFold
  ts.id = allData[[i]]$xtest_IDs
  
  
  ###------------------------------------------------------------------------------------------------------------------------###
  ##  select fold data for train and test splilt ####
  ###------------------------------------------------------------------------------------------------------------------------###
  trainData = td[,-1]
  ytrain = factor(td[,1])
  classes = as.character(unique(ytrain))
  trainIDs = td.id

  
  ###------------------------------------------------------------------------------------------------------------------------###
  ## retrieve test data from kfold partition ####
  ###------------------------------------------------------------------------------------------------------------------------###
  testData = ts[,-1]
  ytest = factor(ts[,1])
  testIDS =ts.id

  
  ###------------------------------------------------------------------------------------------------------------------------###
  ## Apply Data Augmentation ####
  ###------------------------------------------------------------------------------------------------------------------------###
  if(augData){
    
    ## define weights
    tt1 = table(ytrain)
    if(balanceClasses){
      pts.sam = ceiling(max(tt1)/min(tt1))
    }
    
    tt = table(dat[,1])
    ttt = clo(1/tt1)
    tt = data.frame(Status = names(tt),weight = ttt)
    
    ## Aug Data
    ad = simCompData(comp = trainData,labels = ytrain,ptPerSample = pts.sam,fct = 2,alFactor = 1,modthres = modularityThreshold)
    simLabels = data.frame(Status = ad$ClassLabels)
    tt2 = table(simLabels)
    simDat = data.frame(Status = simLabels,ad$simData)
    
    ## BAlance Classes
    if(balanceClasses){
      ## Sample New augmented Points
      ph = data.frame()
      for(c in classes){
        ovrsamp = max(tt1)/tt1[c]
        if(ovrsamp==1){
          nsamp = round(percentMaj * tt1[c])
          ph = rbind(ph,sample_n(simDat[simLabels$Status==c,],size = nsamp))
        }else{
          bool = tt1[c]*pts.sam - max(tt1)
          if(bool==0){
            nsamp=tt1[c]*pts.sam
          }else{
            nsamp = tt1[c]*pts.sam - bool
          }
          ph = rbind(ph,sample_n(simDat[simLabels$Status==c,],size = nsamp))
        }
        
      }
    }else{
      ph = sample_n(simDat,size = numSamples)
    }
   
    

    td = ph[,-1]
    yt = ph[,1]
    table(yt)
    
    ## Redefine Augmented Training Set
    trainData = rbind(trainData,td)
    ytrain = c(as.character(ytrain),as.character(yt))
    table(ytrain)
  }
  

  # ##ROSE
  # # rose_train <- ROSE::ROSE(Status ~ ., data  = data.frame(clr(trainData),Status = ytrain))$data 
  # # trainData = clrInv( rose_train[,-ncol(rose_train)])
  # # ytrain = rose_train[,ncol(rose_train)]
  # 
  
  #S SMOTE
  # rose_train <- DMwR::SMOTE(Status ~ ., data  = data.frame(clr(trainData),Status = ytrain))
  # trainData = clrInv( rose_train[,-ncol(rose_train)])
  # ytrain = rose_train[,ncol(rose_train)]

  ###------------------------------------------------------------------------------------------------------------------------###
  ## Train Model ####
  ###------------------------------------------------------------------------------------------------------------------------###
  nf=5
  mdls = trainML_Models(trainLRs =  data.frame(clr(trainData)),testLRs = data.frame(clr(testData)),
                        ytrain = factor(ytrain),
                        cvMethod = "repeatedCV",mtry_ =  round(sqrt(ncol(trainData[,-1]))),
                        numFolds = nf,numRepeats = 3,
                        y_test = factor(ytrain),testIDs = ts.id,
                        bagModels = F,sampleSize = round(percentBag*nrow(trainData1)),
                        models = c("ranger","pam","svmRadial","glmnet")) 
  
  
  mdls$performance
  models = unique(mdls$predictionMatrix$model)
  
  
  ## get performance AUC
  perf = data.frame()
  for(m in models){
    ph = mdls$predictionMatri %>% 
      filter(model == m )
    
    testLabels = ph[,Outcome]
    probs = ph[,classes]
    rocobj = pROC::multiclass.roc(testLabels,probs)
    
    ll = MLmetrics::MultiLogLoss(y_pred = probs,y_true = testLabels)
    
    perf = rbind(perf,data.frame(Seed = sd,Fold = i,Model = m, 
                                 AUC = as.numeric(pROC::auc(rocobj)),LogLoss = ll))
  }
  
  

  
 
   
  
  ###------------------------------------------------------------------------------------------------------------------------###
  ## Append Output Parm Data  ####
  ###------------------------------------------------------------------------------------------------------------------------###
  sensOptPerformance = rbind(sensOptPerformance,perf)
  View(sensOptPerformance)
  
  
}

allp = sensOptPerformance %>% 
  group_by(Model) %>% 
  summarise_all(.funs = mean)
allp
mean(allp$AUC)
