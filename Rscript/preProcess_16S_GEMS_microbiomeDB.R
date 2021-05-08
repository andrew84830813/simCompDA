###################################################
#### GEMS
####################################################
metaData = data.frame(read_tsv(file = "Data/GEMS.16s_DADA2.sample_details.tsv"))
otuData = data.frame(read_tsv(file = "Data/GEMS.16s_DADA2.taxon_abundance.tsv"))

#define caseIDs
metaData = metaData %>% 
  filter(Case.or.control.subject=="Control")
otuData = subset(otuData,select = c("X1",metaData$X1))

## Sep
otuData1 = separate(otuData,col = 1,into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species","fhfhf"),remove = F,sep = ";")

## retain taxa idenfitfed through family level
otuData1 = otuData1[!is.na(otuData1$Family),]
otuData1 = otuData1[,-9:-2]
otuData1[is.na(otuData1)] = 0
otuData1 = gather(otuData1,"SampleID","Counts",2:ncol(otuData1))
otuData1 = otuData1 %>% 
  group_by(X1,SampleID) %>% 
  summarise(counts = sum(Counts))
otuData1 = spread(otuData1,key = "X1","counts",fill = 0)
xx = str_split(colnames(otuData1[,-1]),pattern = ";",simplify = T)
xx = data.frame(xx) %>% 
  unite("taxa", c(5,6,7), remove = T,sep = "|")
colnames(otuData1) = c("Status",xx$taxa)

write_csv(otuData1,path = "Output/16S_GEMS_microbiomeDB.csv")



###################################################
#### BMI Healthy
####################################################
otu_df = read.delim(file = "Data/crc_baxter.otu_table.100.denovo.rdp_assigned.txt")
sample_metadata = read.delim(file = "Data/crc_baxter.metadata.txt")
sample_metadata$Sample_Name_s = paste0("X",sample_metadata$Sample_Name_s)
diseaseStatus = sample_metadata %>% 
  dplyr::select(Sample_Name_s,DiseaseState) %>% 
  dplyr::rename(X.SampleID = Sample_Name_s)

#parse columns
otu_data = separate(otu_df,1,sep = ";s__",into = c("Taxa"))
otu_data = gather(otu_data,key = "X.SampleID",value = "Counts",2:ncol(otu_data))
otu_data = otu_data %>% 
  group_by(Taxa,X.SampleID) %>% 
  summarise(Counts = sum(Counts))

#Remove poorly designated Taxa
removeDesignations = c("k__;p__;c__;o__;f__;g__",
                       "k__Archaea;p__;c__;o__;f__;g__" ,
                       "k__Bacteria;p__;c__;o__;f__;g__"  )
otu_data = otu_data %>% 
  filter(!Taxa%in%removeDesignations)

# transform sample x feature matrix
otu_data = spread(otu_data,key = "Taxa",value = "Counts",fill = 0) 
otu_data = right_join(diseaseStatus,otu_data)

#get healthy BMI (18.5  - 24.9)
otu_data = otu_data %>% 
  filter(DiseaseState=="H") %>% 
  dplyr::select(-X.SampleID)

colnames(otu_data)[1] = "Status"
df = otu_data
write_csv(otuData1,path = "Output/16S_healthyBMI_microbiomeHD.csv")
##########################################################################################






###################################################
## WGS
##################################################
## VatanenT_2016
load("Data/curatedMetaGenome_VatanenT_2016-control.Rda")
df = test$otuTable
pt.MetaData = test$pt_MetaData
df = df[pt.MetaData$study_condition=="control",]
write_csv(x = df,path = "Output/WGS_controls_vataneen2016.csv")  

## SchirmerM_2016
load("Data/curatedMetaGenome_SchirmerM_2016-control.Rda")
df = test$otuTable
pt.MetaData = test$pt_MetaData
pt.MetaData1 = pt.MetaData %>% 
  filter(smoker=="no")
df = df[pt.MetaData$smoker=="no",]
write_csv(x = df,path = "Output/WGS_controls_SchirmerM_2016.csv")  


load("Data/curatedMetaGenome_HMP_2012-oralCavity.Rda")
df = test$otuTable
pt.MetaData = test$pt_MetaData
write_csv(x = df,path = "Output/WGS_oralCavity_HMP2012.csv")  






##***************************************************
## ml simulations
##***************************************************

###########################################-
#Curated Metagenome (Vincent CDI) #######
###########################################-
load("Data/curatedMetaGenome_Vincent_2016CDI.Rda")
df = test$otuTable
pt.MetaData = test$pt_MetaData
#day0 samples
d0 = pt.MetaData %>% 
  filter(days_from_first_collection==0)
keep = rownames(df)%in%d0$subjectID
df = df[keep,]

df = df %>%
  mutate(Status = if_else(Status==1,"pos","neg")) 
table(df$Status)
write_csv(x = df,path = "Output/WGS_CDI_VINCENT2016.csv")  
###########################################-



#######################################################################################
#load data
otu_df = read.delim(file = "Data/ibd_gevers_2014.otu_table.100.denovo.rdp_assigned")
sample_metadata = read.delim(file = "Data/ibd_gevers_2014.metadata.txt")
diseaseStatus = sample_metadata %>% 
  dplyr::select(sample,DiseaseState) 
colnames(diseaseStatus)[1] = "X.SampleID"

#parse columns
otu_data = separate(otu_df,1,sep = ";s__",into = c("Taxa"))
otu_data = gather(otu_data,key = "X.SampleID",value = "Counts",2:ncol(otu_data))
otu_data = otu_data %>% 
  group_by(Taxa,X.SampleID) %>% 
  summarise(Counts = sum(Counts))

#Remove poorly designated Taxa
removeDesignations = c("k__;p__;c__;o__;f__;g__",
                       "k__Archaea;p__;c__;o__;f__;g__" ,
                       "k__Bacteria;p__;c__;o__;f__;g__"  )
otu_data = otu_data %>% 
  filter(!Taxa%in%removeDesignations)

# transform sample x feature matrix
otu_data = spread(otu_data,key = "Taxa",value = "Counts",fill = 0)
otu_data$X.SampleID = as.character(otu_data$X.SampleID)
diseaseStatus$X.SampleID  =as.character(diseaseStatus$X.SampleID)

otu_data = right_join(diseaseStatus,otu_data)
otu_data = otu_data[,-1]
colnames(otu_data)[1] = "Status"


df = otu_data
colnames(df) = c("Status",paste0("V",1:ncol(df[,-1])))

df = df %>% 
  filter(Status!="NA")
#Assign random labels
table(df$Status)
write_csv(x = df,path = "Output/16S_IBS-unnbalanced_microbiomeHD.csv")
##########################################################################################



############################################
## IHMP Data
############################################
#OTU Data
otu_df = data.frame(read_csv(file = "Data/ihmp16_Data.csv"))
#Sample info
sample_metadata = data.frame(read_csv(file = "Data/ihmp_IBD_metadata.csv"))
sample_metadata$External.ID = paste0("X",sample_metadata$External.ID)
# sample_metadata = sample_metadata %>% 
#   filter(diagnosis=="CD")

#gather data
otu_data = gather(otu_df,key = "External.ID",value = 'Counts',2:ncol(otu_df))
otu_data = otu_data %>% 
  group_by(External.ID,taxonomy) %>% 
  summarise(Counts = sum(Counts)) %>%
  spread(key = "taxonomy",value = "Counts")
bool = as.logical(colSums(otu_data[,-1])!=0)
otu_data = otu_data[,c(T,bool)]

## Join Data
otu_data = right_join(sample_metadata,otu_data)
otu_data = otu_data %>% 
  filter(!is.na(diagnosis))
df = otu_data[,-6:-1]
table(df$diagnosis)
write_csv(x = df,path = "Output/16S_IBDmulticlass_ihmp.csv")


##############################
## bench disease
###########################
fname = "YeZ_2018"
test = getData2(fname = fname,returnCounts = T,site = "stool",
                binaryConditions = NULL,
                positiveClass = NULL)

#Example Data that can be retrieved
pt.MetaData = test$pt_MetaData
df = test$otuTable
df$Status = str_replace_all(df$Status,pattern = "-","")
table(df$Status)
experimentData( test$Info )
abstract( test$Info )
write_csv(x = df,path = "Output/WGS_BenchetDisease_curatedMetagenome.csv")




### feng wgs multiclass
fname = "FengQ_2015"
fname = "ThomasAM_2018a"
fname = "YeZ_2018"
test = getData2(fname = fname,returnCounts = T,site = "stool",
                binaryConditions = NULL,
                positiveClass = NULL)

#Example Data that can be retrieved
pt.MetaData = test$pt_MetaData
df = test$otuTable
df$Status = str_replace_all(df$Status,pattern = "-","")
table(df$Status)
experimentData( test$Info )
abstract( test$Info )

