#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
######################### 09 CHECK RESULTS ##########################
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# Set the working directory
setwd("/imppc/labs/dnalab/share/PRS/09_Check_results/")

# Libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(pROC)
library(caret)
library(e1071)
library(grid)
library(readxl) # read_xlsx
library(openxlsx)


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
##### Check the higher score GCAT participants conditions  005 ######
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

#General score in all individuals
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --pheno-name pheno ",
              " --score  ../08_PRS_advanced/outputs/snps_score_0.05_1e-06_YES", 
              " --out prs_005 "))

#Read score and order 
general_score = fread("prs_005.sscore")
general_score %>% arrange(desc(SCORE1_AVG))

# Density plot to check
density(general_score$SCORE1_AVG)

# Order the scores
general_score = general_score %>% arrange(desc(SCORE1_AVG))

# Filter only the GCAT
general_score = general_score[which(general_score$pheno == 1),]
head(general_score)

# Create the high risk group
# Filter the GCATs with score higher than 0.4 
highrisk = general_score[which(general_score$SCORE1_AVG > 0.42),]
dim(highrisk) # 108 

# Check
length(which(general_score$SCORE1_AVG>0.4))/nrow(general_score)* 100 #107
length(which(general_score$SCORE1_AVG>6))/nrow(general_score)* 100 # 1
length(which(general_score$SCORE1_AVG>4))/nrow(general_score)* 100 # 5
length(which(general_score$SCORE1_AVG>2))/nrow(general_score)* 100 # 11


png("graphs/005_general_GCAT.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
ggplot(general_score,aes((SCORE1_AVG))) + 
  labs(title = "GCAT SCORES in PRS maf 0.05") + xlab("Score") + ylab("") +
  geom_histogram(bins = 100) + geom_vline(xintercept = 0.42) + theme_classic()
dev.off()


png("graphs/005_general_GCAT_highrisk.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
ggplot(highrisk,aes((SCORE1_AVG))) + 
  labs(title = "GCAT highrisk SCORES in PRS maf 0.05") + xlab("Score") + ylab("") +
  geom_histogram(bins = 100) + theme_classic()
dev.off()


png("graphs/005_general_GCAT_highrisk_density.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
ggplot(highrisk,aes((SCORE1_AVG))) + 
  labs(title = "GCAT highrisk SCORES in PRS maf 0.05") + xlab("Score") + ylab("") +
  geom_density() + theme_classic()
dev.off()


# Read gcat
gcat_participants = readRDS("/imppc/labs/dnalab/share/PRS/GCAT_participants/gcat_participants_2019-03-06.rds")
head(gcat_participants)

# Read conditions
conditions = fread("/imppc/labs/dnalab/share/PRS/GCAT_participants/conditions.csv")
head(conditions)

# Merge the conditions and gcat_participants
conditions = left_join(conditions,gcat_participants %>% dplyr::select(entity_id,GENOTYPED_SAMPLE))

# Filter the GCAT participants
table(general_score$pheno)

# Take the higher 50 score participants
risk_score = highrisk

# Take the lower 50 score participants
#risk_score = risk_score %>% arrange(SCORE1_AVG) %>% slice(1:50)

# Filter the information from GCAT of this risk patients
gcat_risk = conditions[conditions$GENOTYPED_SAMPLE %in% risk_score$IID,]

head(gcat_risk)
dim(gcat_risk)
length(unique(gcat_risk$GENOTYPED_SAMPLE)) #108

table(gcat_risk$link)

cancer = gcat_risk[which(gcat_risk$link == "CANCER1"),]

risk_score[which(risk_score$IID %in% cancer$GENOTYPED_SAMPLE),]$SCORE1_AVG
# 0.671134 0.455574 0.420865

# Check the cancer on GCAT
table(conditions$link)
# we have CANCER1, CANCER2 and CANCER3
cancer_conditions = conditions[which(conditions$link == "CANCER1" |
                                       conditions$link == "CANCER2"|
                                       conditions$link == "CANCER3"),]
# We have only the cancer types 
table(cancer_conditions$link)

head(cancer_conditions)

head(general_score)

dim(cancer_conditions) # 322  10
colnames(general_score)[2] = "GENOTYPED_SAMPLE"

cancer_conditions = left_join(cancer_conditions, general_score)
dim(cancer_conditions) # 322  15
head(cancer_conditions)

length(unique(cancer_conditions$GENOTYPED_SAMPLE))#76 individuals

# Witout score and without id
sum(is.na(cancer_conditions$SCORE1_AVG))
sum(is.na(cancer_conditions$GENOTYPED_SAMPLE))

# Order the cancer conditions
cancer_conditions = cancer_conditions %>% arrange(entity_id)

# Delete the ones not genotyped (without genotuped_sample)
cancer_conditions = cancer_conditions[which(!is.na(cancer_conditions$GENOTYPED_SAMPLE)),]
dim(cancer_conditions) # 83 15
length(unique(cancer_conditions$GENOTYPED_SAMPLE))# 75


## Check the prediction for this 75 samples

# Read the score and change column name
general_score_70 = fread("../08_PRS_advanced/outputs/prs_0.05_1e-06_YES_0.07_70.sscore")
colnames(general_score_70)[3] = "general_70"

# Delete na values
general_score = general_score_70 %>% filter(!is.na(general_70))

#boxplot(log10(general_score$SCORE1_AVG)~general_score$general_70)

# Change the 1 and 2 values for 0 and 1 (the ones we need to develop the model)
general_score$general_70 = ifelse(general_score$general_70==1,0,1)

# indicate the characteristics of the train function
ctrl <- trainControl(method = "repeatedcv", number = 5,
                     savePredictions = TRUE)

# Calculate the linear model with the 70% 
lm_general = train(general_70~SCORE1_AVG,data=general_score,method="glm", family="binomial",
                   trControl = ctrl)

# Read the scores file and change the column name 
general_score_30 = fread(paste0("../08_PRS_advanced/outputs/prs_0.05_1e-06_YES_0.07_30.sscore"))
colnames(general_score_30)[3] = "general_30"

# Filter the NAs
general_score_30 = general_score_30 %>% filter(!is.na(general_30))

# To make the linear model (lm) we need to change the 0 and 1 
# because by default the function understand 0 as case and 1 as control 
general_score_30$general_30 = ifelse(general_score_30$general_30==2,1,0)

# Make the prediction with the model, with the 30%
my_pred_general = predict(lm_general,general_score_30,type = "raw")

# We need to put the threshold, in our case we use different thresholds
my_pred_general = ifelse(my_pred_general>0.07,1,0)

# Make a contingency table to visualize the results
my_tab = table(prediction = my_pred_general,
               general = general_score_30$general_30)

my_pred_general = (as.matrix(my_pred_general))
general_score_30 = cbind(general_score_30,my_pred_general)
head(general_score_30)
dim(general_score_30)
table(general_score_30$general_30) # 1473  432 -- > 0 is control (GCAT)
table(general_score_30$general_30, general_score_30$my_pred_general) #33 GCAT classified as case

# filter the control cases 
general_score_30_GCAT = general_score_30[which(general_score_30$general_30 == 0),]
dim(general_score_30_GCAT) # 1473    7

# filter by predicted as 1 (case)
general_score_30_GCAT_cancer = general_score_30_GCAT[which(general_score_30_GCAT$my_pred_general == 1),]
dim(general_score_30_GCAT_cancer) # 33with bad classification

# Now filter the ones with cancer 
length(which(general_score_30_GCAT_cancer$IID %in% cancer_conditions$GENOTYPED_SAMPLE)) # ANY
length(unique(cancer_conditions$GENOTYPED_SAMPLE))

# Cancer in 30%
length(which(general_score_30$IID %in% cancer_conditions$GENOTYPED_SAMPLE)) # 27
# Which ones:
cancer_in_30 = cancer_conditions[which(cancer_conditions$GENOTYPED_SAMPLE %in% general_score_30$IID),]
dim(cancer_in_30) # 1 colon and 10 breast 

cancer_in_30 = cancer_in_30[which(cancer_in_30$desc == "Malignant neoplasm of colon" |
                                    cancer_in_30$desc == "Malignant neoplasm of breast"),]
dim(cancer_in_30) # 1 colon and 10 breast 


# Save the information:
gcat_risk[is.na(gcat_risk)] <- 0 # substitute NAs
write.xlsx(x = gcat_risk, file = "GCAT/cancer_GCAT.xlsx",
           sheetName = "GCAT information", row.names = FALSE)


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
##### Check the higher score GCAT participants conditions  001 ######
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

## Check the prediction for this 75 samples

# Read the score and change column name
general_score_70 = fread("../08_PRS_advanced/outputs/prs_0.01_1e-04_YES_0.05_70.sscore")
colnames(general_score_70)[3] = "general_70"

# Delete na values
general_score = general_score_70 %>% filter(!is.na(general_70))

#boxplot(log10(general_score$SCORE1_AVG)~general_score$general_70)

# Change the 1 and 2 values for 0 and 1 (the ones we need to develop the model)
general_score$general_70 = ifelse(general_score$general_70==1,0,1)

# indicate the characteristics of the train function
ctrl <- trainControl(method = "repeatedcv", number = 5,
                     savePredictions = TRUE)

# Calculate the linear model with the 70% 
lm_general = train(general_70~SCORE1_AVG,data=general_score,method="glm", family="binomial",
                   trControl = ctrl)

# Read the scores file and change the column name 
general_score_30_001 = fread(paste0("../08_PRS_advanced/outputs/prs_0.01_1e-04_YES_0.05_30.sscore"))
colnames(general_score_30_001)[3] = "general_30"

# Filter the NAs
general_score_30_001 = general_score_30_001 %>% filter(!is.na(general_30))

# To make the linear model (lm) we need to change the 0 and 1 
# because by default the function understand 0 as case and 1 as control 
general_score_30_001$general_30 = ifelse(general_score_30_001$general_30==2,1,0)

# Make the prediction with the model, with the 30%
my_pred_general_001 = predict(lm_general,general_score_30_001,type = "raw")

# We need to put the threshold, in our case we use different thresholds
my_pred_general_001 = ifelse(my_pred_general_001>0.05,1,0)

# Make a contingency table to visualize the results
my_tab = table(prediction = my_pred_general_001,
               general = general_score_30_001$general_30)

my_pred_general_001 = (as.matrix(my_pred_general_001))
general_score_30_001 = cbind(general_score_30_001,my_pred_general_001)
head(general_score_30_001)
dim(general_score_30_001)
table(general_score_30_001$general_30) # 1473  432 -- > 0 is control (GCAT)
table(general_score_30_001$general_30, general_score_30_001$my_pred_general_001) 
#60 GCAT classified as case

# filter the control cases 
general_score_30_001 = general_score_30_001[which(general_score_30_001$general_30 == 0),]
dim(general_score_30_001) # 1473    7

# filter by predicted as 1 (case)
general_score_30_GCAT_cancer_001 = general_score_30_001[which(general_score_30_001$my_pred_general_001 == 1),]
dim(general_score_30_GCAT_cancer_001) # 60 with bad classification

# Now filter the ones with cancer 
length(which(general_score_30_GCAT_cancer_001$IID %in% cancer_conditions$GENOTYPED_SAMPLE)) # ANY

# Cancer in 30%
length(which(general_score_30_001$IID %in% cancer_conditions$GENOTYPED_SAMPLE)) # 27
# Which ones:
cancer_in_30_001 = cancer_conditions[which(cancer_conditions$GENOTYPED_SAMPLE %in% general_score_30_001$IID),]
dim(cancer_in_30_001) # 1 colon and 10 breast 
table(cancer_in_30_001$desc)

cancer_in_30_001 = cancer_in_30_001[which(cancer_in_30_001$desc == "Malignant neoplasm of colon" |
                                            cancer_in_30_001$desc == "Malignant neoplasm of breast"),]
dim(cancer_in_30_001) # 1 colon and 10 breast 

### Same results as with 005 


# Save the information:
gcat_risk[is.na(gcat_risk)] <- 0 # substitute NAs
write.xlsx(x = gcat_risk, file = "GCAT/cancer_GCAT.xlsx",
           sheetName = "GCAT information", row.names = FALSE)


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
##### Check the higher score GCAT participants conditions  0001 ######
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

## Check the prediction for this 75 samples

# Read the score and change column name
general_score_70 = fread("../08_PRS_advanced/outputs/prs_0.001_0.05_YES_0.05_70.sscore")
colnames(general_score_70)[3] = "general_70"

# Delete na values
general_score = general_score_70 %>% filter(!is.na(general_70))

#boxplot(log10(general_score$SCORE1_AVG)~general_score$general_70)

# Change the 1 and 2 values for 0 and 1 (the ones we need to develop the model)
general_score$general_70 = ifelse(general_score$general_70==1,0,1)

# indicate the characteristics of the train function
ctrl <- trainControl(method = "repeatedcv", number = 5,
                     savePredictions = TRUE)

# Calculate the linear model with the 70% 
lm_general = train(general_70~SCORE1_AVG,data=general_score,method="glm", family="binomial",
                   trControl = ctrl)

# Read the scores file and change the column name 
general_score_30_0001 = fread(paste0("../08_PRS_advanced/outputs/prs_0.001_0.05_YES_0.05_30.sscore"))
colnames(general_score_30_0001)[3] = "general_30"

# Filter the NAs
general_score_30_0001 = general_score_30_0001 %>% filter(!is.na(general_30))

# To make the linear model (lm) we need to change the 0 and 1 
# because by default the function understand 0 as case and 1 as control 
general_score_30_0001$general_30 = ifelse(general_score_30_0001$general_30==2,1,0)

# Make the prediction with the model, with the 30%
my_pred_general_0001 = predict(lm_general,general_score_30_0001,type = "raw")

# We need to put the threshold, in our case we use different thresholds
my_pred_general_0001 = ifelse(my_pred_general_0001>0.05,1,0)

# Make a contingency table to visualize the results
my_tab = table(prediction = my_pred_general_0001,
               general = general_score_30_0001$general_30)

my_pred_general_0001 = (as.matrix(my_pred_general_0001))
general_score_30_0001 = cbind(general_score_30_0001,my_pred_general_0001)
head(general_score_30_0001)
dim(general_score_30_0001)
table(general_score_30_0001$general_30) # 1473  432 -- > 0 is control (GCAT)
table(general_score_30_0001$general_30, general_score_30_0001$my_pred_general_0001) 
#88 GCAT classified as case

# filter the control cases 
general_score_30_0001 = general_score_30_0001[which(general_score_30_0001$general_30 == 0),]
dim(general_score_30_0001) # 1473    7

# filter by predicted as 1 (case)
general_score_30_GCAT_cancer_0001 = general_score_30_0001[which(general_score_30_0001$my_pred_general_0001 == 1),]
dim(general_score_30_GCAT_cancer_0001) # 88 with bad classification

# Now filter the ones with cancer 
length(which(general_score_30_GCAT_cancer_0001$IID %in% cancer_conditions$GENOTYPED_SAMPLE)) # ANY

# Cancer in 30%
length(which(general_score_30_0001$IID %in% cancer_conditions$GENOTYPED_SAMPLE)) # 27
# Which ones:
cancer_in_30_0001 = cancer_conditions[which(cancer_conditions$GENOTYPED_SAMPLE %in% general_score_30_0001$IID),]
dim(cancer_in_30_0001) # 1 colon and 10 breast 
table(cancer_in_30_0001$desc)

cancer_in_30_0001 = cancer_in_30_0001[which(cancer_in_30_0001$desc == "Malignant neoplasm of colon" |
                                            cancer_in_30_0001$desc == "Malignant neoplasm of breast"),]
dim(cancer_in_30_0001) # 1 colon and 9 breast 

### Same results as with 005 

ids = cancer_in_30_0001$entity_id

conditions_cancer = conditions[which(conditions$entity_id %in% ids),]
dim(conditions_cancer)
length(unique(conditions_cancer$entity_id)) # 10 individuals
conditions_cancer = conditions_cancer %>% arrange(conditions_cancer$entity_id)

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
#################  CREATE A TABLE WITH CANCER INFO ##################
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# Read gcat
gcat_participants = readRDS("/imppc/labs/dnalab/share/PRS/GCAT_participants/gcat_participants_2019-03-06.rds")
colnames(gcat_participants) 
# we are interested in: entity_id, GENOTYPED_SAMPLE, BIRTH_YEAR, Gender,DISEASE

# Read conditions
conditions = fread("/imppc/labs/dnalab/share/PRS/GCAT_participants/conditions.csv")
colnames(conditions)

# Merge the conditions and gcat_participants
conditions = left_join(conditions,gcat_participants %>% dplyr::select(entity_id, GENOTYPED_SAMPLE, AGE, BIRTH_YEAR, GENDER))
head(conditions)

# Filter the cancer people
conditions_allcancer_ids = conditions[which(conditions$link == "CANCER1" |
                                          conditions$link == "CANCER2"|
                                          conditions$link == "CANCER3"),]
dim(conditions_allcancer_ids) # 322
conditions_allcancer_ids = unique(conditions_allcancer_ids$entity_id)
length(conditions_allcancer_ids) # 293 individuals 

conditions_allcancer = conditions[which(conditions$entity_id %in% conditions_allcancer_ids),]
dim(conditions_allcancer)
head(conditions_allcancer)

conditions_allcancer = conditions_allcancer %>% arrange(entity_id)
conditions_allcancer$dat = NULL


# Now mark with the 30%
group = rep(0,nrow(conditions_allcancer))
group = ifelse(conditions_allcancer$entity_id %in% ids, 30, 70)
conditions_allcancer$group = group

# Check 
test = conditions_allcancer[which(conditions_allcancer$entity_id %in% ids),]
table(test$group)

# Filter the genotyped
conditions_allcancer = conditions_allcancer[which(!is.na(conditions_allcancer$GENOTYPED_SAMPLE)),]

# check 
length(unique(conditions_allcancer$entity_id))#75

#ADD age
conditions_allcancer$age_in_2019 = 2019-as.integer(conditions_allcancer$BIRTH_YEAR)

# age-diagnostic
age_diagnostic = conditions_allcancer$dat0
age_diagnostic[is.na(age_diagnostic)] <- 0 # substitute NAs
age_diagnostic = tstrsplit(age_diagnostic, "-")[1]
age_diagnostic = age_diagnostic[[1]]
age_diagnostic = conditions_allcancer$age_in_2019-(2019-as.integer(age_diagnostic))
age_diagnostic[age_diagnostic < 0] <- 0

conditions_allcancer$age_diagnosis = age_diagnostic

# ARRANGE INFO
conditions_allcancer = conditions_allcancer %>% dplyr::select(entity_id,GENOTYPED_SAMPLE,group,
                                                              AGE,age_in_2019, GENDER, 
                                                              age_diagnosis,link,  
                                                              desc, antecedents,relation)

head(conditions_allcancer)
colnames(conditions_allcancer) = c("Entity_id", "Genotyped_id", "Group", "Age_start_GCAT", 
                                   "Patient_Age", "Gender", "Age_diagnosis", "Disease",
                                   "Disease_description", "Antecedents", "Relation")

# Save the information:
conditions_allcancer[is.na(conditions_allcancer)] <- 0 # substitute NAs

write.xlsx(x = conditions_allcancer, file = "GCAT/GCAT_allcancerinformation.xlsx",
           sheetName = "GCAT information", row.names = FALSE)


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
##########  CREATE A TABLE WITH GCAT CLASIFIED AS CASE INFO #########
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# First we need to identify the bad classified individuals of each model 

head(general_score_30_GCAT)

# We need to filter only the bad classified
table(general_score_30_GCAT$my_pred_general)# 33
bad_class_005 = general_score_30_GCAT[which(general_score_30_GCAT$my_pred_general == 1),]$IID

# Same for 001 model
head(general_score_30_001)
dim(general_score_30_001)
table(general_score_30_001$general_30)
table(general_score_30_001$my_pred_general_001)# 60 
bad_class_001 = general_score_30_001[which(general_score_30_001$my_pred_general_001 == 1),]$IID

# Same for 0001 model
head(general_score_30_0001)
dim(general_score_30_0001)
table(general_score_30_0001$general_30)
table(general_score_30_0001$my_pred_general_0001)# 88
bad_class_0001 = general_score_30_0001[which(general_score_30_0001$my_pred_general_0001 == 1),]$IID

# Check the IIDs
length(which(bad_class_005 %in% bad_class_0001)) #27/33
length(which(bad_class_001 %in% bad_class_0001)) #47/60

# Merge all the IDDs
bad_allids = c(bad_class_005,bad_class_001, bad_class_0001)
length(bad_allids)

bad_allids = unique(bad_allids)
length(bad_allids) # 101

# Read gcat
gcat_participants = readRDS("/imppc/labs/dnalab/share/PRS/GCAT_participants/gcat_participants_2019-03-06.rds")
colnames(gcat_participants) 
# we are interested in: entity_id, GENOTYPED_SAMPLE, BIRTH_YEAR, Gender,DISEASE

# Read conditions
conditions = fread("/imppc/labs/dnalab/share/PRS/GCAT_participants/conditions.csv")
colnames(conditions)

# Merge the conditions and gcat_participants
conditions = left_join(conditions,gcat_participants %>% dplyr::select(entity_id, GENOTYPED_SAMPLE, AGE, BIRTH_YEAR, GENDER))
head(conditions)

conditions_bad_class = conditions[which(conditions$GENOTYPED_SAMPLE %in% bad_allids),]
dim(conditions_bad_class)
length(unique(conditions_bad_class$entity_id)) #101

conditions_bad_class = conditions_bad_class %>% arrange(entity_id)

#ADD age
conditions_bad_class$age_in_2019 = 2019-as.integer(conditions_bad_class$BIRTH_YEAR)

# age-diagnostic
age_diagnostic = conditions_bad_class$dat0
age_diagnostic[is.na(age_diagnostic)] <- 0 # substitute NAs
age_diagnostic = tstrsplit(age_diagnostic, "-")[1]
age_diagnostic = age_diagnostic[[1]]
age_diagnostic = conditions_bad_class$age_in_2019-(2019-as.integer(age_diagnostic))
age_diagnostic[age_diagnostic < 0] <- 0

conditions_bad_class$age_diagnosis = age_diagnostic

# ARRANGE INFO
conditions_bad_class = conditions_bad_class %>% dplyr::select(entity_id,GENOTYPED_SAMPLE,
                                                              AGE,age_in_2019, GENDER, 
                                                              age_diagnosis,link,  
                                                              desc, antecedents,relation)

head(conditions_bad_class)
colnames(conditions_bad_class) = c("Entity_id", "Genotyped_id", "Age_start_GCAT", 
                                   "Patient_Age", "Gender", "Age_diagnosis", "Disease",
                                   "Disease_description", "Antecedents", "Relation")

# Create columns
conditions_bad_class$model005 = rep("Control",nrow(conditions_bad_class))
conditions_bad_class$model001 = rep("Control",nrow(conditions_bad_class))
conditions_bad_class$model0001 = rep("Control",nrow(conditions_bad_class))

# Add the models columns
for (r in 1:nrow(conditions_bad_class)){
  id = conditions_bad_class$Genotyped_id[r]
  if (id %in% bad_class_005){
    conditions_bad_class$model005[r] = "Case"
  }
  if (id %in% bad_class_001){
    conditions_bad_class$model001[r] = "Case"
  }
  if (id %in% bad_class_0001){
    conditions_bad_class$model0001[r] = "Case"
  }
}

head(conditions_bad_class)

# Check
length(bad_class_005) # 33
length(bad_class_001) #60
length(bad_class_0001) #88

test = conditions_bad_class[!duplicated(conditions_bad_class$Entity_id),]
table(test$model005)
table(test$model001)
table(test$model0001)

# Save the information:
conditions_bad_class[is.na(conditions_bad_class)] <- 0 # substitute NAs

write.xlsx(x = conditions_bad_class, file = "GCAT/GCAT_badclassifiedinformation.xlsx",
           sheetName = "GCAT information", row.names = FALSE)


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
####################  CHECK CANCERS CASES ON 30% ####################
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# Read the information
conditions_allcancer = read_xlsx("GCAT/GCAT_allcancerinformation.xlsx")
conditions_allcancer = as.data.frame(conditions_allcancer)
head(conditions_allcancer)

# Filter the ones in the 30
str(conditions_allcancer)
table(conditions_allcancer$Group)
conditions_allcancer = conditions_allcancer %>% filter(Group == 30)
dim(conditions_allcancer) #115
length(unique(conditions_allcancer$Entity_id))# 10

# Check it
ids = unique(conditions_allcancer$Entity_id)
diseases = c("Malignant neoplasm of colon", 
             "Malignant neoplasm of breast",
             "Malignant neoplasm of bladder",
             "Malignant neoplasm of prostate")
ids_with_antecedents = c()
for (i in ids){
  patient = conditions_allcancer %>% 
    filter(conditions_allcancer$Entity_id == i & conditions_allcancer$Antecedents== 0)
  family = conditions_allcancer %>% 
    filter(conditions_allcancer$Entity_id == i & conditions_allcancer$Antecedents== 1)
  if (length(which(family$Disease_description %in% diseases)) > 0){
    ids_with_antecedents = c(ids_with_antecedents, i)
  }
}

GCAT_withantecedents= conditions_allcancer %>% filter(conditions_allcancer$Entity_id %in% ids_with_antecedents)

length(unique(GCAT_withantecedents$Entity_id)) # 5 patients with antecedents
#check
GCAT_withantecedents[,8:11]
conditions_allcancer %>% filter(conditions_allcancer$Entity_id =="=E00251532000721")

write.xlsx(x = GCAT_withantecedents, file = "GCAT/GCAT_cancer30%_withantecedents.xlsx",
           sheetName = "GCAT information", row.names = FALSE)

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
##############  CHECK ANTECEDENTS IN THE BADCLASSIFIED ##############
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
conditions_bad_class = read_xlsx("GCAT/GCAT_badclassifiedinformation.xlsx")
conditions_bad_class = as.data.frame(conditions_bad_class)
head(conditions_bad_class)

ids = unique(conditions_bad_class$Entity_id) #101
diseases = unique(conditions_bad_class[grep("Malignant neoplasm", conditions_bad_class$Disease_description),]$Disease_description)
ids_with_antecedents = c()
for (i in ids){
  patient = conditions_bad_class %>% filter(conditions_bad_class$Entity_id == i)
  if (length(which(patient$Disease_description %in% diseases)) > 0){
    ids_with_antecedents = c(ids_with_antecedents, i)
  }
}


GCAT_withantecedents_badclass= conditions_bad_class %>% 
  filter(conditions_bad_class$Entity_id %in% ids_with_antecedents)

length(unique(GCAT_withantecedents_badclass$Entity_id)) # 48 patients with antecedents
# from 101 patients 48 have antecedents

GCAT_withantecedents_badclass  = GCAT_withantecedents_badclass %>% arrange(Entity_id)

GCAT_withantecedents_badclass$color = rep(0,nrow(GCAT_withantecedents_badclass))
GCAT_withantecedents_badclass$cancer = rep(0,nrow(GCAT_withantecedents_badclass))

prev = "NO"
count = 0
for (i in 1:nrow(GCAT_withantecedents_badclass)){
  print(i)
  actual = GCAT_withantecedents_badclass$Entity_id[i]
  if (GCAT_withantecedents_badclass$Disease_description[i] %in% diseases){
    GCAT_withantecedents_badclass$cancer[i] = "CANCER"
    print("YES")
  }
  if (actual != prev){
    count = count + 1
  }
  prev = actual
  GCAT_withantecedents_badclass$color[i] = count
}
table(GCAT_withantecedents_badclass$cancer)
dim(GCAT_withantecedents_badclass)
GCAT_withantecedents_badclass = GCAT_withantecedents_badclass[which(GCAT_withantecedents_badclass$Disease_description != 0),]

write.xlsx(x = GCAT_withantecedents_badclass, file = "GCAT/GCAT_badclassified30%_withantecedents.xlsx",
           sheetName = "GCAT information", row.names = FALSE)

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
##############  CHECK ANTECEDENTS IN THE BADCLASSIFIED ##############
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# Read gcat
gcat_participants = readRDS("/imppc/labs/dnalab/share/PRS/GCAT_participants/gcat_participants_2019-03-06.rds")
colnames(gcat_participants) 
# we are interested in: entity_id, GENOTYPED_SAMPLE, BIRTH_YEAR, Gender,DISEASE
dim(gcat_participants)
table(gcat_participants$SNP_IMPUTED)
gcat_participants = gcat_participants[which(gcat_participants$SNP_IMPUTED == T),]
dim(gcat_participants)
# Read conditions
conditions = fread("/imppc/labs/dnalab/share/PRS/GCAT_participants/conditions.csv")
colnames(conditions)

# Merge the conditions and gcat_participants
conditions = left_join(conditions,gcat_participants %>% dplyr::select(entity_id, GENOTYPED_SAMPLE, AGE, BIRTH_YEAR, GENDER))
head(conditions)

# Take only the genotyped
conditions = conditions[which(!is.na(conditions$GENOTYPED_SAMPLE)),]
dim(conditions) #47968
length(unique(conditions$entity_id)) #4988

ids = unique(conditions_bad_class$Entity_id) #101

# filter the ones not in ids
conditions_wellclas = conditions[-which(conditions$entity_id %in% ids),]
dim(conditions_wellclas) # 46972 
length(unique(conditions_wellclas$entity_id)) # 4887

# but we want only the ones in the 30%
general_score_30 = fread(paste0("../08_PRS_advanced/outputs/prs_0.05_1e-06_YES_0.07_30.sscore"))
colnames(general_score_30)[3] = "general_30"
head(general_score_30)

general_score_30 = general_score_30[-which(is.na(general_score_30$general_30)),]
dim(general_score_30)# 1905
table(general_score_30$general_30)
general_score_30 = general_score_30 %>% filter(general_30 == 1)
dim(general_score_30)# 1473

# Filter the ones in the 30
conditions_wellclas = conditions_wellclas[which(conditions_wellclas$GENOTYPED_SAMPLE %in% general_score_30$IID),]
dim(conditions_wellclas) # 13208
length(unique(conditions_wellclas$entity_id))# 1372
     
# order by ID
conditions_wellclas = conditions_wellclas %>% arrange(entity_id)

#ADD age
conditions_wellclas$age_in_2019 = 2019-as.integer(conditions_wellclas$BIRTH_YEAR)

# age-diagnostic
age_diagnostic = conditions_wellclas$dat0
age_diagnostic[is.na(age_diagnostic)] <- 0 # substitute NAs
age_diagnostic = tstrsplit(age_diagnostic, "-")[1]
age_diagnostic = age_diagnostic[[1]]
age_diagnostic = conditions_wellclas$age_in_2019-(2019-as.integer(age_diagnostic))
age_diagnostic[age_diagnostic < 0] <- 0

conditions_wellclas$age_diagnosis = age_diagnostic

# ARRANGE INFO
conditions_wellclas = conditions_wellclas %>% dplyr::select(entity_id,GENOTYPED_SAMPLE,
                                                              AGE,age_in_2019, GENDER, 
                                                              age_diagnosis,link,  
                                                              desc, antecedents,relation)

head(conditions_wellclas)
colnames(conditions_wellclas) = c("Entity_id", "Genotyped_id", "Age_start_GCAT", 
                                   "Patient_Age", "Gender", "Age_diagnosis", "Disease",
                                   "Disease_description", "Antecedents", "Relation")




ids = unique(conditions_wellclas$Entity_id) #101
diseases = unique(conditions_wellclas[grep("Malignant neoplasm", conditions_wellclas$Disease_description),]$Disease_description)
ids_with_antecedents = c()
for (i in ids){
  patient = conditions_wellclas %>% filter(conditions_wellclas$Entity_id == i)
  if (length(which(patient$Disease_description %in% diseases)) > 0){
    ids_with_antecedents = c(ids_with_antecedents, i)
  }
}

length(ids_with_antecedents)# 552
GCAT_withantecedents_wellclass= conditions_wellclas %>% 
  filter(conditions_wellclas$Entity_id %in% ids_with_antecedents)

GCAT_withantecedents_wellclass  = GCAT_withantecedents_wellclass %>% arrange(Entity_id)

GCAT_withantecedents_wellclass$color = rep(0,nrow(GCAT_withantecedents_wellclass))
GCAT_withantecedents_wellclass$cancer = rep(0,nrow(GCAT_withantecedents_wellclass))

prev = "NO"
count = 0
for (i in 1:nrow(GCAT_withantecedents_wellclass)){
  actual = GCAT_withantecedents_wellclass$Entity_id[i]
  if (GCAT_withantecedents_wellclass$Disease_description[i] %in% diseases){
    GCAT_withantecedents_wellclass$cancer[i] = "CANCER"
  }
  if (actual != prev){
    count = count + 1
  }
  prev = actual
  GCAT_withantecedents_wellclass$color[i] = count
}
table(GCAT_withantecedents_wellclass$cancer)# 2929
dim(GCAT_withantecedents_wellclass)
GCAT_withantecedents_wellclass = GCAT_withantecedents_wellclass[which(GCAT_withantecedents_wellclass$Disease_description != 0),]

# Save the information:
GCAT_withantecedents_wellclass[is.na(GCAT_withantecedents_wellclass)] <- 0 # substitute NAs

length(unique(GCAT_withantecedents_wellclass$Entity_id)) 

write.xlsx(x = GCAT_withantecedents_wellclass, file = "GCAT/GCAT_wellclassifiedinformation_fulltable.xlsx",
           sheetName = "GCAT information", row.names = FALSE)

table(GCAT_withantecedents_wellclass$cancer)

GCAT_withantecedents_wellclass = GCAT_withantecedents_wellclass[which(GCAT_withantecedents_wellclass$cancer == "CANCER"),]
dim(GCAT_withantecedents_wellclass) # 820
length(unique(GCAT_withantecedents_wellclass$Entity_id)) # 552

write.xlsx(x = GCAT_withantecedents_wellclass, file = "GCAT/GCAT_wellclassifiedinformation_allcancers.xlsx",
           sheetName = "GCAT information", row.names = FALSE)

#filter the cancers
our_cancers = c("Malignant neoplasm of breast",
                "Malignant neoplasm of colon",
                "Malignant neoplasm of ovary")

GCAT_withantecedents_wellclass = GCAT_withantecedents_wellclass %>% 
  filter(Disease_description %in% our_cancers)

dim(GCAT_withantecedents_wellclass) # 234
length(unique(GCAT_withantecedents_wellclass$Entity_id)) # 200

write.xlsx(x = GCAT_withantecedents_wellclass, file = "GCAT/GCAT_wellclassifiedinformation_ourcancers.xlsx",
           sheetName = "GCAT information", row.names = FALSE)

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
#######################  CHECK CANCERS IN 30% ######################
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# Read gcat
gcat_participants = readRDS("/imppc/labs/dnalab/share/PRS/GCAT_participants/gcat_participants_2019-03-06.rds")
colnames(gcat_participants) 
# we are interested in: entity_id, GENOTYPED_SAMPLE, BIRTH_YEAR, Gender,DISEASE
dim(gcat_participants)
table(gcat_participants$SNP_IMPUTED)
gcat_participants = gcat_participants[which(gcat_participants$SNP_IMPUTED == T),]
dim(gcat_participants)
# Read conditions
conditions = fread("/imppc/labs/dnalab/share/PRS/GCAT_participants/conditions.csv")
colnames(conditions)

# Read catalogue
catalogue = fread("/imppc/labs/dnalab/share/PRS/GCAT_participants/catalogue.csv")
colnames(catalogue)
head(catalogue)
table(catalogue$domain)
catalogue = catalogue %>% filter(catalogue$domain == "Conditions") # filter for conditions
dim(catalogue)
diseases = unique(catalogue[grep("C", catalogue$cod),]$cod) # filter all the cancers codes
diseases = c(diseases, c("199", "224", "365")) # add some cancers not in the C category
diseases = catalogue[which(catalogue$cod %in% diseases),]$desc

# Merge the conditions and gcat_participants
conditions = left_join(conditions,gcat_participants %>% dplyr::select(entity_id, GENOTYPED_SAMPLE, AGE, BIRTH_YEAR, GENDER))
head(conditions)

# Take only the genotyped
conditions = conditions[which(!is.na(conditions$GENOTYPED_SAMPLE)),]
dim(conditions) #47968
length(unique(conditions$entity_id)) #4988

# but we want only the ones in the 30%
general_score_30 = fread(paste0("../08_PRS_advanced/outputs/prs_0.05_1e-06_YES_0.07_30.sscore"))
colnames(general_score_30)[3] = "general_30"
head(general_score_30)

general_score_30 = general_score_30[-which(is.na(general_score_30$general_30)),]
dim(general_score_30)# 1905
table(general_score_30$general_30)
general_score_30 = general_score_30 %>% filter(general_30 == 1)
dim(general_score_30)# 1473 ALL GCATs

# Filter the ones in the 30
conditions_30 = conditions[which(conditions$GENOTYPED_SAMPLE %in% general_score_30$IID),]
dim(conditions_30) # 14204
length(unique(conditions_30$entity_id))# 1473

# order by ID
conditions_30 = conditions_30 %>% arrange(entity_id)

#ADD age
conditions_30$age_in_2019 = 2019-as.integer(conditions_30$BIRTH_YEAR)

# age-diagnostic
age_diagnostic = conditions_30$dat0
age_diagnostic[is.na(age_diagnostic)] <- 0 # substitute NAs
age_diagnostic = tstrsplit(age_diagnostic, "-")[1]
age_diagnostic = age_diagnostic[[1]]
age_diagnostic = conditions_30$age_in_2019-(2019-as.integer(age_diagnostic))
age_diagnostic[age_diagnostic < 0] <- 0

conditions_30$age_diagnosis = age_diagnostic

# ARRANGE INFO
conditions_30 = conditions_30 %>% dplyr::select(entity_id,GENOTYPED_SAMPLE,
                                                            AGE,age_in_2019, GENDER, 
                                                            age_diagnosis,link,  
                                                            desc, antecedents,relation)

head(conditions_30)
colnames(conditions_30) = c("Entity_id", "Genotyped_id", "Age_start_GCAT", 
                                  "Patient_Age", "Gender", "Age_diagnosis", "Disease",
                                  "Disease_description", "Antecedents", "Relation")




ids = unique(conditions_30$Entity_id) #101
ids_with_antecedents = c()
for (i in ids){
  patient = conditions_30 %>% filter(conditions_30$Entity_id == i)
  if (length(which(patient$Disease_description %in% diseases)) > 0){# all that have 1 or more diseases descriptions
    ids_with_antecedents = c(ids_with_antecedents, i)
  }
}

length(ids_with_antecedents)# 552
GCAT_withantecedents_30= conditions_30 %>% 
  filter(conditions_30$Entity_id %in% ids_with_antecedents)

GCAT_withantecedents_30  = GCAT_withantecedents_30 %>% arrange(Entity_id)

GCAT_withantecedents_30$color = rep(0,nrow(GCAT_withantecedents_30))
GCAT_withantecedents_30$cancer = rep(0,nrow(GCAT_withantecedents_30))

prev = "NO"
count = 0
for (i in 1:nrow(GCAT_withantecedents_30)){
  actual = GCAT_withantecedents_30$Entity_id[i]
  if (GCAT_withantecedents_30$Disease_description[i] %in% diseases){
    GCAT_withantecedents_30$cancer[i] = "CANCER"
  }
  if (actual != prev){
    count = count + 1
  }
  prev = actual
  GCAT_withantecedents_30$color[i] = count
}
table(GCAT_withantecedents_30$cancer)# 1036
dim(GCAT_withantecedents_30)
GCAT_withantecedents_30 = GCAT_withantecedents_30[which(GCAT_withantecedents_30$Disease_description != 0),]

# Save the information:
GCAT_withantecedents_30[is.na(GCAT_withantecedents_30)] <- 0 # substitute NAs

length(unique(GCAT_withantecedents_30$Entity_id)) #666

write.xlsx(x = GCAT_withantecedents_30, file = "GCAT/GCAT_30%_allinformation.xlsx",
           sheetName = "GCAT information", row.names = FALSE)

table(GCAT_withantecedents_30$cancer)

GCAT_withantecedents_30 = GCAT_withantecedents_30[which(GCAT_withantecedents_30$cancer == "CANCER"),]
dim(GCAT_withantecedents_30) # 1036
length(unique(GCAT_withantecedents_30$Entity_id)) # 666

write.xlsx(x = GCAT_withantecedents_30, file = "GCAT/GCAT_30%_cancers.xlsx",
           sheetName = "GCAT information", row.names = FALSE)

#filter the cancers
our_cancers = c("Malignant neoplasm of breast",
                "Malignant neoplasm of colon",
                "Malignant neoplasm of ovary")

GCAT_withantecedents_30 = GCAT_withantecedents_30 %>% 
  filter(Disease_description %in% our_cancers)

dim(GCAT_withantecedents_30) # 262
length(unique(GCAT_withantecedents_30$Entity_id)) # 221

write.xlsx(x = GCAT_withantecedents_30, file = "GCAT/GCAT_30%_ourcancers.xlsx",
           sheetName = "GCAT information", row.names = FALSE)


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
#######################  CHECK CANCERS IN 70% #######################
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

length(unique(conditions$entity_id)) #4988

# but we want only the ones in the 70%
general_score_70 = fread(paste0("../08_PRS_advanced/outputs/prs_0.05_1e-06_YES_0.07_70.sscore"))
colnames(general_score_70)[3] = "general_70"
head(general_score_70)

general_score_70 = general_score_70[-which(is.na(general_score_70$general_70)),]
dim(general_score_70)# 4458
table(general_score_70$general_70)
general_score_70 = general_score_70 %>% filter(general_70 == 1)
dim(general_score_70)# 3515 ALL GCATs

# Filter the ones in the 70
conditions_70 = conditions[which(conditions$GENOTYPED_SAMPLE %in% general_score_70$IID),]
dim(conditions_70) # 33764
length(unique(conditions_70$entity_id))# 3515

# order by ID
conditions_70 = conditions_70 %>% arrange(entity_id)

#ADD age
conditions_70$age_in_2019 = 2019-as.integer(conditions_70$BIRTH_YEAR)

# age-diagnostic
age_diagnostic = conditions_70$dat0
age_diagnostic[is.na(age_diagnostic)] <- 0 # substitute NAs
age_diagnostic = tstrsplit(age_diagnostic, "-")[1]
age_diagnostic = age_diagnostic[[1]]
age_diagnostic = conditions_70$age_in_2019-(2019-as.integer(age_diagnostic))
age_diagnostic[age_diagnostic < 0] <- 0

conditions_70$age_diagnosis = age_diagnostic

# ARRANGE INFO
conditions_70 = conditions_70 %>% dplyr::select(entity_id,GENOTYPED_SAMPLE,
                                                AGE,age_in_2019, GENDER, 
                                                age_diagnosis,link,  
                                                desc, antecedents,relation)

head(conditions_70)
colnames(conditions_70) = c("Entity_id", "Genotyped_id", "Age_start_GCAT", 
                            "Patient_Age", "Gender", "Age_diagnosis", "Disease",
                            "Disease_description", "Antecedents", "Relation")




ids = unique(conditions_70$Entity_id) 
ids_with_antecedents = c()
for (i in ids){
  patient = conditions_70 %>% filter(conditions_70$Entity_id == i)
  if (length(which(patient$Disease_description %in% diseases)) > 0){# all that have 1 or more diseases descriptions
    ids_with_antecedents = c(ids_with_antecedents, i)
  }
}

length(ids_with_antecedents)# 1574
GCAT_withantecedents_70= conditions_70 %>% 
  filter(conditions_70$Entity_id %in% ids_with_antecedents)

GCAT_withantecedents_70  = GCAT_withantecedents_70 %>% arrange(Entity_id)

GCAT_withantecedents_70$color = rep(0,nrow(GCAT_withantecedents_70))
GCAT_withantecedents_70$cancer = rep(0,nrow(GCAT_withantecedents_70))

prev = "NO"
count = 0
for (i in 1:nrow(GCAT_withantecedents_70)){
  actual = GCAT_withantecedents_70$Entity_id[i]
  if (GCAT_withantecedents_70$Disease_description[i] %in% diseases){
    GCAT_withantecedents_70$cancer[i] = "CANCER"
  }
  if (actual != prev){
    count = count + 1
  }
  prev = actual
  GCAT_withantecedents_70$color[i] = count
}
table(GCAT_withantecedents_70$cancer)# 2406
dim(GCAT_withantecedents_70)
GCAT_withantecedents_70 = GCAT_withantecedents_70[which(GCAT_withantecedents_70$Disease_description != 0),]

# Save the information:
GCAT_withantecedents_70[is.na(GCAT_withantecedents_70)] <- 0 # substitute NAs

length(unique(GCAT_withantecedents_70$Entity_id)) #1574

write.xlsx(x = GCAT_withantecedents_70, file = "GCAT/GCAT_70%_allinformation.xlsx",
           sheetName = "GCAT information", row.names = FALSE)

table(GCAT_withantecedents_70$cancer)

GCAT_withantecedents_70 = GCAT_withantecedents_70[which(GCAT_withantecedents_70$cancer == "CANCER"),]
dim(GCAT_withantecedents_70) # 2406
length(unique(GCAT_withantecedents_70$Entity_id)) # 1574

write.xlsx(x = GCAT_withantecedents_70, file = "GCAT/GCAT_70%_cancers.xlsx",
           sheetName = "GCAT information", row.names = FALSE)

#filter the cancers
GCAT_withantecedents_70 = GCAT_withantecedents_70 %>% 
  filter(Disease_description %in% our_cancers)

dim(GCAT_withantecedents_70) # 678
length(unique(GCAT_withantecedents_70$Entity_id)) # 582

write.xlsx(x = GCAT_withantecedents_70, file = "GCAT/GCAT_70%_ourcancers.xlsx",
           sheetName = "GCAT information", row.names = FALSE)


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
####### Check that GCAT info from GCAT and Pandora is the same ######
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# Now we want to check that the information of the GCAT participants from Pandora and GCAT is the same
# Read pandora information
load("/imppc/labs/dnalab/share/PRS/02_Pandora_data_preparation/data/AllSamples.RData")
colnames(all.samples)[13] = "varsion"
# filter GCAT in pandora
all.samples = all.samples %>% filter(tissue=="DNA")
# filter the ones we have in GCAT 
gcat_in_pandora = gcat_participants[which(gcat_participants$`PI-2017-01` %in% all.samples$SampleName),]
head(gcat_in_pandora)

head(all.samples)

colnames(gcat_in_pandora)[45] = "SampleName"

# take the samples information
all.samples = left_join(all.samples,gcat_in_pandora %>% dplyr::select(SampleName,GENOTYPED_SAMPLE))

dim(all.samples)

# read the gcat information
psam = fread("../05_GWAS/outputs_plink/first_final_subset.psam")

all.samples = all.samples[all.samples$GENOTYPED_SAMPLE %in% psam$IID,]

concordance = NULL
# Check "duplicateds" (presents in pandora and GCAT)
for(i in 1:177){
  
  psam[psam$IID %in% c(all.samples$sampleId[i],all.samples$GENOTYPED_SAMPLE[i]),1:2] %>% 
    write.table("duplicated.txt",row.names = F,col.names = F,quote = F)
  
  
  system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
                " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
                " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
                " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
                " --keep duplicated.txt ", 
                " --recode vcf ", 
                " --out duplicated_check"))
  
  my_vcf = fread("duplicated_check.vcf")
  
  head(my_vcf)
  
  # Take both samples from each individual
  sample1 = as.character(t(my_vcf[,10]))
  sample2 = as.character(t(my_vcf[,11]))
  
  aux = sum(sample1==sample2)/nrow(my_vcf)
  
  # Check if they are similar
  concordance = c(concordance,aux)
  
}

# Check the results with a histogram and summary
hist(concordance)
summary(concordance)
# the results are very similar

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
################### Check the maf 0.05 PRS model  ###################
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# We want to check the SNPs included in this PRS model
# Read the GWAS
gwas = fread("../08_PRS_advanced/outputs/gwas.pheno_70.glm.firth")
dim(gwas) # 13023 variants

# Apply the filters 
# filter for the LD variants in ld_in
ld_in = fread("../08_PRS_advanced/outputs/LD.prune.in",header = F)
gwas = gwas[gwas$ID %in% ld_in$V1,] 
dim(gwas) # 8689 variants
# filter by mad and pvalue
gwas = gwas %>% filter(A1_CASE_FREQ>0.05 & P<1e-06)
dim(gwas) # 41 variants

# These 41 variants are the final variants included on the model 
head(gwas)

# We filter to keep only the variables we are interested in
variants_inmodel = gwas[,1:5]
head(variants_inmodel)

# Create the chr pos
variants_inmodel$chr_pos = paste0(variants_inmodel$`#CHROM`, "_", variants_inmodel$POS)
head(variants_inmodel)

# Now we want to know the gene per each variant
genes_information = fread("../04_Merge_datasets/output_data/panels_subset.txt")
# Create the chr pos
genes_information = genes_information[,1:2]
head(genes_information)
dim(genes_information) # 13747 variants information 


# Check if selected variants in genes_information
sum(genes_information$chr_pos %in% variants_inmodel$chr_pos) # 41

# filter the genes_information we are interested in
genes_information = genes_information[which(genes_information$chr_pos %in% variants_inmodel$chr_pos),]

# Now merge to obtain all the information 
variants_inmodel = merge(variants_inmodel,genes_information)
variants_inmodel  = variants_inmodel %>% dplyr::select(chr_pos,`#CHROM`,POS,ID, Gene)

head(variants_inmodel)

length(unique(variants_inmodel$Gene)) # 34 genes in the model 

genes_names = unique(variants_inmodel$Gene)
genes_names = genes_names[order(genes_names)]
genes_names

# Read the studied genes information 
studied_genes = read_xlsx("../Common_variants_colonbreastoverlap literature2018.xlsx", sheet = "genes")
head(studied_genes)
dim(studied_genes) #308  19

# First we want to know if this genes are in GWAs PRS (other studies report this genes)
# filter from studied_genes only the ones we are interested in
studied_genes = studied_genes[which(studied_genes$`Genes common variants GWA analysis` %in% variants_inmodel$Gene),]
head(studied_genes)
dim(studied_genes) #35
# Change columns names
colnames(studied_genes)[6] = "Genes"
colnames(studied_genes)[3] = "GWAS_PRS"
colnames(studied_genes)[4] = "Neoplasia_overall"

#delete not interesting columns
studied_genes = studied_genes[,-2]
studied_genes = studied_genes[,-4]

# As data frame
studied_genes = as.data.frame(studied_genes)
head(studied_genes)

length(unique(studied_genes$Genes)) # We have one duplication, delete it
studied_genes = studied_genes[-6,]
dim(studied_genes)

studied_genes[is.na(studied_genes)] <- 0
studied_genes[which(studied_genes$`MAMA_v2_20-02-19` == "MAMA_v2_20-02-19"),]$`MAMA_v2_20-02-19` <- "MAMA"
studied_genes[which(studied_genes$`OVARIO_v1_11-11-16` == "OVARIO_v1_11-11-16"),]$`OVARIO_v1_11-11-16` <- "OVARIO"
studied_genes[which(studied_genes$`CMOH_v1_11-11-16` == "CMOH_v1_11-11-16"),]$`CMOH_v1_11-11-16` <- "OVARIO Y MAMA"
studied_genes[which(studied_genes$`POLIPOSIS_v1_11-11-16` == "POLIPOSIS_v1_11-11-16"),]$`POLIPOSIS_v1_11-11-16` <- "COLON"
studied_genes[which(studied_genes$`CCHNP_v1_11-11-16` == "CCHNP_v1_11-11-16"),]$`CCHNP_v1_11-11-16` <- "COLON"
studied_genes[which(studied_genes$`COLON_JOVEN_v1_11-11-16` == "COLON_JOVEN_v1_11-11-16"),]$`COLON_JOVEN_v1_11-11-16` <- "COLON"
studied_genes[which(studied_genes$`LI-FRAUMENI_v1_11-11-16` == "LI-FRAUMENI_v1_11-11-16"),]$`LI-FRAUMENI_v1_11-11-16` <- "LI-FRAUMENI"
studied_genes[which(studied_genes$`MELANOMA_v3_17-05-18` == "MELANOMA_v3_17-05-18"),]$`MELANOMA_v3_17-05-18` <- "MELANOMA"
studied_genes[which(studied_genes$`GASTRICO_v1_11-11-16` == "GASTRICO_v1_11-11-16"),]$`GASTRICO_v1_11-11-16` <- "GASTRICO"
studied_genes[which(studied_genes$`Erdhein Chester_v1_28-03-19` == "Erdhein Chester_v1_28-03-19"),]$`Erdhein Chester_v1_28-03-19` <- "SANGRE"
studied_genes[which(studied_genes$`RENAL_v1_11-11-16` == "RENAL_v1_11-11-16"),]$`RENAL_v1_11-11-16` <- "RENAL"
studied_genes[which(studied_genes$`PANCREAS_v1_9-9-17` == "PANCREAS_v1_9-9-17"),]$`PANCREAS_v1_9-9-17` <- "PANCREAS"
studied_genes[which(studied_genes$`PROSTATA_v2_9-9-17` == "PROSTATA_v2_9-9-17"),]$`PROSTATA_v2_9-9-17` <- "PROSTATA"

# Save the information:
write.xlsx(x = studied_genes, file = "genes_inPRS_005.xlsx",
           sheetName = "Genes information", row.names = FALSE)

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
################### Check the maf 0.01 PRS model  ###################
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# We want to check the SNPs included in this PRS model
# Read the GWAS
gwas = fread("../08_PRS_advanced/outputs/gwas.pheno_70.glm.firth")
dim(gwas) # 13023 variants

# Apply the filters 
# filter for the LD variants in ld_in
ld_in = fread("../08_PRS_advanced/outputs/LD.prune.in",header = F)
gwas = gwas[gwas$ID %in% ld_in$V1,] 
dim(gwas) # 8689 variants
# filter by mad and pvalue
gwas = gwas %>% filter(A1_CASE_FREQ>0.01 & P<0.0001)
dim(gwas) # 105 variants

# These 105 variants are the final variants included on the model 
head(gwas)

# We filter to keep only the variables we are interested in
variants_inmodel = gwas[,1:5]
head(variants_inmodel)

# Create the chr pos
variants_inmodel$chr_pos = paste0(variants_inmodel$`#CHROM`, "_", variants_inmodel$POS)
head(variants_inmodel)

# Now we want to know the gene per each variant
genes_information = fread("../04_Merge_datasets/output_data/panels_subset.txt")
# Create the chr pos
genes_information = genes_information[,1:2]
head(genes_information)
dim(genes_information) # 13747 variants information 


# Check if selected variants in genes_information
sum(genes_information$chr_pos %in% variants_inmodel$chr_pos) # 105

# filter the genes_information we are interested in
genes_information = genes_information[which(genes_information$chr_pos %in% variants_inmodel$chr_pos),]

# Now merge to obtain all the information 
variants_inmodel = merge(variants_inmodel,genes_information)
variants_inmodel  = variants_inmodel %>% dplyr::select(chr_pos,`#CHROM`,POS,ID, Gene)

head(variants_inmodel)

length(unique(variants_inmodel$Gene)) # 66 genes in the model 

genes_names = unique(variants_inmodel$Gene)
genes_names = genes_names[order(genes_names)]
genes_names

# Read the studied genes information 
studied_genes = read_xlsx("../Common_variants_colonbreastoverlap literature2018.xlsx", sheet = "genes")
head(studied_genes)
dim(studied_genes) #308  19

# First we want to know if this genes are in GWAs PRS (other studies report this genes)
# filter from studied_genes only the ones we are interested in
studied_genes = studied_genes[which(studied_genes$`Genes common variants GWA analysis` %in% variants_inmodel$Gene),]
head(studied_genes)
dim(studied_genes) # 67
# Change columns names
colnames(studied_genes)[6] = "Genes"
colnames(studied_genes)[3] = "GWAS_PRS"
colnames(studied_genes)[4] = "Neoplasia_overall"

#delete not interesting columns
studied_genes = studied_genes[,-2]
studied_genes = studied_genes[,-4]

# As data frame
studied_genes = as.data.frame(studied_genes)
head(studied_genes)

length(unique(studied_genes$Genes)) # We have one duplication, delete it
studied_genes = studied_genes[-11,]
dim(studied_genes)

studied_genes[is.na(studied_genes)] <- 0

studied_genes[is.na(studied_genes)] <- 0
studied_genes[which(studied_genes$`MAMA_v2_20-02-19` == "MAMA_v2_20-02-19"),]$`MAMA_v2_20-02-19` <- "MAMA"
studied_genes[which(studied_genes$`OVARIO_v1_11-11-16` == "OVARIO_v1_11-11-16"),]$`OVARIO_v1_11-11-16` <- "OVARIO"
studied_genes[which(studied_genes$`CMOH_v1_11-11-16` == "CMOH_v1_11-11-16"),]$`CMOH_v1_11-11-16` <- "OVARIO Y MAMA"
studied_genes[which(studied_genes$`POLIPOSIS_v1_11-11-16` == "POLIPOSIS_v1_11-11-16"),]$`POLIPOSIS_v1_11-11-16` <- "COLON"
studied_genes[which(studied_genes$`CCHNP_v1_11-11-16` == "CCHNP_v1_11-11-16"),]$`CCHNP_v1_11-11-16` <- "COLON"
studied_genes[which(studied_genes$`COLON_JOVEN_v1_11-11-16` == "COLON_JOVEN_v1_11-11-16"),]$`COLON_JOVEN_v1_11-11-16` <- "COLON"
studied_genes[which(studied_genes$`LI-FRAUMENI_v1_11-11-16` == "LI-FRAUMENI_v1_11-11-16"),]$`LI-FRAUMENI_v1_11-11-16` <- "LI-FRAUMENI"
studied_genes[which(studied_genes$`MELANOMA_v3_17-05-18` == "MELANOMA_v3_17-05-18"),]$`MELANOMA_v3_17-05-18` <- "MELANOMA"
studied_genes[which(studied_genes$`GASTRICO_v1_11-11-16` == "GASTRICO_v1_11-11-16"),]$`GASTRICO_v1_11-11-16` <- "GASTRICO"
studied_genes[which(studied_genes$`Erdhein Chester_v1_28-03-19` == "Erdhein Chester_v1_28-03-19"),]$`Erdhein Chester_v1_28-03-19` <- "SANGRE"
studied_genes[which(studied_genes$`RENAL_v1_11-11-16` == "RENAL_v1_11-11-16"),]$`RENAL_v1_11-11-16` <- "RENAL"
studied_genes[which(studied_genes$`PANCREAS_v1_9-9-17` == "PANCREAS_v1_9-9-17"),]$`PANCREAS_v1_9-9-17` <- "PANCREAS"
studied_genes[which(studied_genes$`PROSTATA_v2_9-9-17` == "PROSTATA_v2_9-9-17"),]$`PROSTATA_v2_9-9-17` <- "PROSTATA"

# Save the information:
write.xlsx(x = studied_genes, file = "genes_inPRS_001.xlsx",
           sheetName = "Genes information", row.names = FALSE)


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
################### Check the maf 0.001 PRS model  ##################
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# We want to check the SNPs included in this PRS model
# Read the GWAS
gwas = fread("../08_PRS_advanced/outputs/gwas.pheno_70.glm.firth")
dim(gwas) # 13023 variants

# Apply the filters 
# filter for the LD variants in ld_in
ld_in = fread("../08_PRS_advanced/outputs/LD.prune.in",header = F)
gwas = gwas[gwas$ID %in% ld_in$V1,] 
dim(gwas) # 8689 variants
# filter by mad and pvalue
gwas = gwas %>% filter(A1_CASE_FREQ>0.001 & P<0.05)
dim(gwas) # 844 variants

# These 844 variants are the final variants included on the model 
head(gwas)

# We filter to keep only the variables we are interested in
variants_inmodel = gwas[,1:5]
head(variants_inmodel)

# Create the chr pos
variants_inmodel$chr_pos = paste0(variants_inmodel$`#CHROM`, "_", variants_inmodel$POS)
head(variants_inmodel)

# Now we want to know the gene per each variant
genes_information = fread("../04_Merge_datasets/output_data/panels_subset.txt")
# Create the chr pos
genes_information = genes_information[,1:2]
head(genes_information)
dim(genes_information) # 13747 variants information 


# Check if selected variants in genes_information
sum(genes_information$chr_pos %in% variants_inmodel$chr_pos) # 859

# filter the genes_information we are interested in
genes_information = genes_information[which(genes_information$chr_pos %in% variants_inmodel$chr_pos),]

# Now merge to obtain all the information 
variants_inmodel = merge(variants_inmodel,genes_information)
variants_inmodel  = variants_inmodel %>% dplyr::select(chr_pos,`#CHROM`,POS,ID, Gene)

head(variants_inmodel)

length(unique(variants_inmodel$Gene)) # 129 genes in the model 

genes_names = unique(variants_inmodel$Gene)
genes_names = genes_names[order(genes_names)]
genes_names

# Read the studied genes information 
studied_genes = read_xlsx("../Common_variants_colonbreastoverlap literature2018.xlsx", sheet = "genes")
head(studied_genes)
dim(studied_genes) #308  19

# First we want to know if this genes are in GWAs PRS (other studies report this genes)
# filter from studied_genes only the ones we are interested in
studied_genes = studied_genes[which(studied_genes$`Genes common variants GWA analysis` %in% variants_inmodel$Gene),]
head(studied_genes)
dim(studied_genes) # 129
# Change columns names
colnames(studied_genes)[6] = "Genes"
colnames(studied_genes)[3] = "GWAS_PRS"
colnames(studied_genes)[4] = "Neoplasia_overall"

#delete not interesting columns
studied_genes = studied_genes[,-2]
studied_genes = studied_genes[,-4]

# As data frame
studied_genes = as.data.frame(studied_genes)
head(studied_genes)

length(unique(studied_genes$Genes)) # We have one duplication, delete it
studied_genes = studied_genes[-18,]
dim(studied_genes)

studied_genes[is.na(studied_genes)] <- 0

studied_genes[is.na(studied_genes)] <- 0
studied_genes[which(studied_genes$`MAMA_v2_20-02-19` == "MAMA_v2_20-02-19"),]$`MAMA_v2_20-02-19` <- "MAMA"
studied_genes[which(studied_genes$`OVARIO_v1_11-11-16` == "OVARIO_v1_11-11-16"),]$`OVARIO_v1_11-11-16` <- "OVARIO"
studied_genes[which(studied_genes$`CMOH_v1_11-11-16` == "CMOH_v1_11-11-16"),]$`CMOH_v1_11-11-16` <- "OVARIO Y MAMA"
studied_genes[which(studied_genes$`POLIPOSIS_v1_11-11-16` == "POLIPOSIS_v1_11-11-16"),]$`POLIPOSIS_v1_11-11-16` <- "COLON"
studied_genes[which(studied_genes$`CCHNP_v1_11-11-16` == "CCHNP_v1_11-11-16"),]$`CCHNP_v1_11-11-16` <- "COLON"
studied_genes[which(studied_genes$`COLON_JOVEN_v1_11-11-16` == "COLON_JOVEN_v1_11-11-16"),]$`COLON_JOVEN_v1_11-11-16` <- "COLON"
studied_genes[which(studied_genes$`LI-FRAUMENI_v1_11-11-16` == "LI-FRAUMENI_v1_11-11-16"),]$`LI-FRAUMENI_v1_11-11-16` <- "LI-FRAUMENI"
studied_genes[which(studied_genes$`MELANOMA_v3_17-05-18` == "MELANOMA_v3_17-05-18"),]$`MELANOMA_v3_17-05-18` <- "MELANOMA"
studied_genes[which(studied_genes$`GASTRICO_v1_11-11-16` == "GASTRICO_v1_11-11-16"),]$`GASTRICO_v1_11-11-16` <- "GASTRICO"
studied_genes[which(studied_genes$`Erdhein Chester_v1_28-03-19` == "Erdhein Chester_v1_28-03-19"),]$`Erdhein Chester_v1_28-03-19` <- "SANGRE"
studied_genes[which(studied_genes$`RENAL_v1_11-11-16` == "RENAL_v1_11-11-16"),]$`RENAL_v1_11-11-16` <- "RENAL"
studied_genes[which(studied_genes$`PANCREAS_v1_9-9-17` == "PANCREAS_v1_9-9-17"),]$`PANCREAS_v1_9-9-17` <- "PANCREAS"
studied_genes[which(studied_genes$`PROSTATA_v2_9-9-17` == "PROSTATA_v2_9-9-17"),]$`PROSTATA_v2_9-9-17` <- "PROSTATA"

# Save the information:
write.xlsx(x = studied_genes, file = "genes_inPRS_0001.xlsx",
           sheetName = "Genes information", row.names = FALSE)


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
################# Create the final table by GENE  ###################
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# Read the gene information from the three models 
model005 = as.data.frame(read_xlsx("genes_inPRS_005.xlsx"))
model001 = as.data.frame(read_xlsx("genes_inPRS_001.xlsx"))
model0001 = as.data.frame(read_xlsx("genes_inPRS_0001.xlsx"))

# Check 
head(model005)
head(model001)
head(model0001)
dim(model005) # 34 genes
dim(model001) # 66 genes
dim(model0001) # 128 genes

# Check the genes
sum(model005$Genes %in% model0001$Genes) # all 
sum(model001$Genes %in% model0001$Genes) # all

#order the genes
model005 = model005[order(model005$Genes),]
model001 = model001[order(model001$Genes),]
model0001 = model0001[order(model0001$Genes),]

m001 = rep(0,128)
m001 = ifelse((allgenes_models$Genes %in% model001$Genes), "YES", "NO")
m005 = rep(0,128)
m005 = ifelse((allgenes_models$Genes %in% model005$Genes), "YES", "NO")

# create the new dataframe
allgenes_models = data.frame(model0001$Genes,
                             rep(0,128),  # chr
                             rep("YES",128),
                             m001,#model001
                             m005)#model005
                            
#change colnames 
colnames(allgenes_models) = c("Genes","chr", "model0001", 
                              "model001", "model005")                            
head(allgenes_models)


# Check
allgenes_models[1:10,1:5]
sum(allgenes_models$model0001 == "YES") #128
sum(allgenes_models$model001 == "YES") # 66
sum(allgenes_models$model005 == "YES") # 34

# Check the chromosome
head(genes_information)
head(allgenes_models)
for (i in 1:nrow(allgenes_models)){
  chr = as.character(tstrsplit(genes_information[genes_information$Gene == allgenes_models$Genes[i],2][1],"_")[1])
  allgenes_models$chr[i] = chr
}

# take the information to add 
table_info = studied_genes %>% dplyr::select(Genes,`MAMA_v2_20-02-19`,`OVARIO_v1_11-11-16`,
                                             `CMOH_v1_11-11-16`,`POLIPOSIS_v1_11-11-16`,
                                             `CCHNP_v1_11-11-16`,`COLON_JOVEN_v1_11-11-16`,
                                             `LI-FRAUMENI_v1_11-11-16`,`MELANOMA_v3_17-05-18`,
                                             `GASTRICO_v1_11-11-16`,`Erdhein Chester_v1_28-03-19`,
                                             `RENAL_v1_11-11-16`,`PANCREAS_v1_9-9-17`,
                                             `PROSTATA_v2_9-9-17`, AHG, GWAS_PRS,
                                             Neoplasia_overall)

head(table_info)
colnames(table_info)

head(table_info[,14:16])

for (i in 1:nrow(table_info)){
  if (table_info$GWAS_PRS[i] == "GWAs"){
    table_info$GWAS_PRS[i] = table_info$Neoplasia_overall[i]
  }
}

table(table_info$GWAS_PRS)
table_info$Neoplasia_overall = NULL
# change colnames
colnames(table_info) = c("Genes","Breast","Ovarian","CMOH: Ovarian and Breast", "Poliposis: Colon",
                         "CCHNP: Colon", "Young Colon", "Li-Fraumeni", "Melanoma", "Gastric",
                         "Erdhein Chester: Blood", "Renal", "Pancreas", "Prostate", "AHG",
                         "GWAS_PRS")
head(table_info)

# Change the values of the cancer types columns
for (j in 2:14){
  table_info[,j] = ifelse(table_info[,j] == 0, 0, 1)
}

# We have to merge allgenes_models and table_info
head(allgenes_models)
head(table_info)
allgenes_models$Genes = as.character(allgenes_models$Genes)
allgenes_models_final = left_join(allgenes_models,table_info)
head(allgenes_models_final)

# Save the information:
write.xlsx(x = allgenes_models_final, file = "genes_inPRS_allmodels.xlsx",
           sheetName = "Genes information", row.names = FALSE)

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
################## Create the final table by SNP ####################
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

head(allgenes_models)
head(table_info)
head(gwas)
head(genes_information)
dim(genes_information) # 1 SNP missing 

gwas = fread("/imppc/labs/dnalab/share/PRS/08_PRS_advanced/outputs/gwas.pheno_70.glm.firth")
dim(gwas) # 13023 variants
# filter for the LD variants in ld_in
ld_in = fread("../08_PRS_advanced/outputs/LD.prune.in",header = F)
gwas = gwas[gwas$ID %in% ld_in$V1,] # 8689 variants

# We have to create the gwas subset for each model
gwas_0001 = gwas %>% filter(A1_CASE_FREQ>0.001 & P<0.05) # 844 variants
gwas_001 = gwas %>% filter(A1_CASE_FREQ>0.01 & P<0.0001) # 105 variants
gwas_005 = gwas %>% filter(A1_CASE_FREQ>0.05 & P<1e-06) # 41 variants

#check

gwas_0001[which(!gwas_0001$chr_pos %in% genes_information$chr_pos),]

#844 unique IDs
length(unique(gwas_0001[!is.na(gwas_0001$ID),]$ID))

# Start with gene name and chr_pos
# Create the chromosome and position for each variable
genes_information$Chromosome = tstrsplit(genes_information$chr_pos,"_")[1]
genes_information$Position = tstrsplit(genes_information$chr_pos,"_")[2]
# check
head(genes_information)
# Need to create the chr_pos for gwas
gwas$chr_pos = paste0(gwas$`#CHROM`,"_",gwas$POS)

# Create the CI variable

# Now select the interesting information 
gwas = gwas %>% dplyr::select(chr_pos, ID, P,OR,L95,U95,A1_CASE_FREQ,A1_CTRL_FREQ)
gwas_0001 = gwas_0001 %>% dplyr::select(ID, P,OR,L95,U95,A1_CASE_FREQ,A1_CTRL_FREQ)
gwas_001 = gwas_001 %>% dplyr::select(ID, P,OR,L95,U95,A1_CASE_FREQ,A1_CTRL_FREQ)
gwas_005 = gwas_005 %>% dplyr::select(ID, P,OR,L95,U95,A1_CASE_FREQ,A1_CTRL_FREQ)

# Merge the data
final_table = left_join(genes_information,gwas)
head(final_table)
colnames(final_table)

# Create the CI variable
final_table$CI = paste0(final_table$L95,"-",final_table$U95)
final_table$L95 = NULL
final_table$U95 = NULL
colnames(final_table)

# Now add the table_info
colnames(table_info)[1] = "Gene"
final_table= left_join(final_table,table_info, by= "Gene")
head(final_table)

# Now we have to add the models
colnames(allgenes_models)[1] = "Gene"
allgenes_models$chr = NULL
final_table = left_join(final_table,allgenes_models)
final_table$chr_pos = NULL
head(final_table)

# Check the real SNPs for each model 
final_table$model0001 = ifelse(final_table$ID %in% gwas_0001$ID, "Yes", "No")
final_table$model001 = ifelse(final_table$ID %in% gwas_001$ID, "Yes", "No")
final_table$model005 = ifelse(final_table$ID %in% gwas_005$ID, "Yes", "No")


final_table$AHG = ifelse(final_table$AHG == 0, 0, "Yes")

# check 
length(final_table[which(!final_table$ID %in% gwas_0001$ID),]$ID)
length(gwas_0001[which(!gwas_0001$ID %in% final_table$ID),]$ID)

# One missing and repeatitions
length(unique(final_table$ID))
missing = gwas_0001[which(!gwas_0001$ID %in% final_table$ID),] # the missing one
missing$CI = paste0(missing$L95, "-", missing$U95)
missing$L95 = NULL
missing$U95 = NULL
final_table[1,]
# Need: gene, chr, pos, breast, ovarian, (...), model0001, model001 and model005
gwas_missing = gwas[which(gwas$ID == missing$ID),]
gwas_missing$Gene = "PRSS1"
gwas_missing$`#CHROM` = tstrsplit(gwas_missing$chr_pos,"_")[1]
gwas_missing$POS = tstrsplit(gwas_missing$chr_pos,"_")[2]

# take gene chr and pos information
missing_info = gwas_missing %>% dplyr::select(Gene, `#CHROM`, POS)
missing_info = cbind(missing_info, missing) # add 
# take cascades information
missing = table_info[which(table_info$Gene == missing_info$Gene),]
missing$Gene = NULL
missing_info = cbind(missing_info, missing)
# Take models information
missing = allgenes_models[which(allgenes_models$Gene == missing_info$Gene),]
missing$Gene = NULL
missing$model0001 = "Yes"
missing$model001 = "Yes"
missing$model005 = "Yes"
missing_info = cbind(missing_info, missing)
colnames(missing_info) = colnames(final_table)
final_table  = rbind(final_table, missing_info)
dim(final_table)
final_table = final_table %>% arrange((Gene)) # order by gene name

# Delete duplicated
length(unique(final_table$ID)) # 844
final_table = final_table[!duplicated(final_table$ID),]
dim(final_table) # 844  27

# Check the models snps
table(final_table$model0001)
table(final_table$model001)
table(final_table$model005)

# Save the information:
write.xlsx(x = final_table, file = "SNPs_inPRS_allmodels.xlsx",
           sheetName = "Genes information", row.names = FALSE)

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
################## Create the final table by GENE ###################
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

head(final_table)

# Number of genes:
length(unique(final_table$Gene)) #129 genes
genes = unique(final_table$Gene)

final_table_genes = data.frame(matrix(ncol = 27, nrow = 0))
colnames(final_table_genes) = colnames(final_table)
for (i in genes){
  subset = final_table[final_table$Gene == i,]
  aux = c()
  for (j in 1:nrow(subset)){
    aux = c(aux, sum(subset[j,25:27] == "Yes"))  
  }
  final_table_genes = rbind(final_table_genes, subset[which.max(aux),])
}


final_table_genes[is.na(final_table_genes)] <- 0


# Save the information:
final_table_genes[is.na(final_table_genes)] <- 0 # substitute NAs
write.xlsx(x = final_table_genes, file = "genes_inPRS_allmodels.xlsx",
           sheetName = "Genes information", row.names = FALSE)


