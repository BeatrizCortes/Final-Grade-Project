######################################################################################
# 07 PRS                                                                             #
######################################################################################

# Set the working directory
setwd("/imppc/labs/dnalab/share/PRS/07_PRS")

# Libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(pROC)

# Read the data 
psam = fread("../05_GWAS/outputs_plink/first_final_subset.psam")

head(psam)

# number of control and case by group
table(psam$pheno_70) # 3515  943 
table(psam$ovarian_70) # 3506  295
table(psam$breast_70) # 3492  584 
table(psam$colon_70) # 3489  219  

# Load the general variants 
general_variants2 = fread("../06_GWAS/70_general_variants.txt")
ovarian_variants2 = fread("../06_GWAS/70_ovarian_variants.txt")
breast_variants2 = fread("../06_GWAS/70_breast_variants.txt")
colon_variants2 = fread("../06_GWAS/70_colon_variants.txt")

# Check the number of correlated SNPs for each group
dim(general_variants2) # 385
dim(ovarian_variants2) # 169
dim(breast_variants2) # 297
dim(colon_variants2) # 124

# Create the directories for the results
dir.create("General")
dir.create("Ovarian")
dir.create("Breast")
dir.create("Colon")


##############################################################
# GENERAL ANALYSIS                                           #
##############################################################

# Create the file with the beast coefficients, contain: ID, Allele and freq
general_variants2[,c(3,6,9)] %>% write.table("General/general_coefficients.txt",
                                            row.names = F,col.names = F,quote = F)

# Calculate the general score with PLINK with a maf of 0.01
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --pheno-name pheno_70 ",
              " --score  General/general_coefficients.txt ", 
              " --maf 0.01 ", 
              " --out General/general_score"))

# Read the output .sscore file 
general_score = fread("General/general_score.sscore")
head(general_score)
colnames(general_score)[3] = "general_70"
table(general_score$general_70) # 3515  943 
sum(is.na(general_score$general_70)) # 2111

# Filter the NAs
general_score = general_score %>% filter(!is.na(general_70))

# Check the results of both groups 1 = control, 2 = case 
tapply(general_score$SCORE1_AVG,INDEX = general_score$general_70,summary)

# Plot: boxplot 
boxplot(log10(general_score$SCORE1_AVG)~general_score$general_70)

# Kruskal test 
kruskal.test(general_score$SCORE1_AVG~general_score$general_70)

# Plot: boxplot
png("General/general_boxplot_70.png")
ggplot(general_score,aes(as.factor(general_70),SCORE1_AVG)) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: boxplot
png("General/general_boxplotlog10_70.png")
ggplot(general_score,aes(as.factor(general_70),log10(SCORE1_AVG))) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: histogram
png("General/general_histogram_70.png")
ggplot(general_score,aes(SCORE1_AVG, fill= as.factor(general_70))) + 
  geom_histogram(alpha = 0.5)+
  theme_classic()
dev.off()

# Plot: histogram 
png("General/general_density_70.png")
ggplot(general_score,aes(log10(SCORE1_AVG), fill= as.factor(general_70))) + 
  geom_density(alpha = 0.5)+
  labs(x = "Score", y = "Density log10", fill= "Groups")+
  scale_fill_discrete(name="Group",
                      breaks=c("0", "1"),
                      labels=c("General Cancer", "Health"))+
  theme_classic()
dev.off()

# As we can see almost all the cases are nead the 0 of score and the cases have higher values 

# To make the linear model (lm) we need to change the 0 and 1 
# because by default the function understand 0 as case and 1 as control 
general_score$general_70 = ifelse(general_score$general_70==2,0,1)

table(general_score$general_70) #   943 3515

# Calculate the linear model with the 70% 
lm = glm(general_70~SCORE1_AVG,data=general_score,family = "binomial")
summary(lm)

# Make the prediction with the model, with the 70%
my_pred = predict(lm,general_score,type = "response")

# We need to put the threshold, in our case 0.95
my_pred = ifelse(my_pred>0.95,1,0)

# Make a contingency table to visualize the results
my_tab = table(prediction = my_pred,
               general = general_score$general_70)
my_tab
#           general
#prediction    0    1
#         0  855 1096
#         1   88 2419

# We have:
# 855: CORRECT cancer as cancer --> 90.66808
# 2419: CORRECT control as health --> 68.81935
# 88: INCORRECT cancer as health --> 9.33192
# 1096: INCORRECT health as cancer  --> 31.18065

# Calculate sensitivity and specificity 
# Sensitivity: percentage of sick people classified as sick
sensitivity = my_tab[1,1]/table(general_score$general_70)[1] # 0.9066808 

# Specificity: percentage of health people classified as health
specificity = my_tab[2,2]/(my_tab[1,2]+my_tab[2,2]) # 0.6881935

# Calculate AUC
auc(my_pred,general_score$general_70)# Area under the curve: 0.7016

# Try to check by thresholding the data
aux = general_score %>% filter(SCORE1_AVG>0.60330) # group 0 + 1st Qu.

table(aux$general_70)


########################################
# Now make the prediction with the 30% #
########################################
head(general_variants)
head(general_score)

# Calculate the general score with PLINK with a maf of 0.01
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --pheno-name pheno_30 ",
              " --score  General/general_coefficients.txt ", 
              " --maf 0.01 ", 
              " --out General/general_score_30"))

# Read the output .sscore file 
general_score_30 = fread("General/general_score_30.sscore")
head(general_score_30)
colnames(general_score_30)[3] = "general_30"
table(general_score_30$general_30) # 1473  432 
sum(is.na(general_score_30$general_30)) #  4664

# Filter the NAs
general_score_30 = general_score_30 %>% filter(!is.na(general_30))

# Check the results of both groups 1 = control, 2 = case 
tapply(general_score_30$SCORE1_AVG,INDEX = general_score_30$general_30,summary)

# Plot: boxplot 
boxplot(log10(general_score_30$SCORE1_AVG)~general_score_30$general_30)

# Kruskal test 
kruskal.test(general_score_30$SCORE1_AVG~general_score_30$general_30)

# Plot: boxplot
png("General/general_boxplot_30.png")
ggplot(general_score_30,aes(as.factor(general_30),SCORE1_AVG)) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: boxplot
png("General/general_boxplotlog10_30.png")
ggplot(general_score_30,aes(as.factor(general_30),log10(SCORE1_AVG))) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: histogram
png("General/general_histogram_30.png")
ggplot(general_score_30,aes(SCORE1_AVG, fill= as.factor(general_30))) + 
  geom_histogram(alpha = 0.5)+
  theme_classic()
dev.off()

# Plot: histogram 
png("General/general_density_30.png")
ggplot(general_score_30,aes(log10(SCORE1_AVG), fill= as.factor(general_30))) + 
  geom_density(alpha = 0.5)+
  labs(x = "Score", y = "Density log10", fill= "Groups")+
  scale_fill_discrete(name="Group",
                      breaks=c("0", "1"),
                      labels=c("general Cancer", "Health"))+
  theme_classic()
dev.off()

# To make the linear model (lm) we need to change the 0 and 1 
# because by default the function understand 0 as case and 1 as control 
general_score_30$general_30 = ifelse(general_score_30$general_30==2,0,1)

table(general_score_30$general_30) #  432 1473 

# Make the prediction with the model, with the 30%
my_pred = predict(lm,general_score_30,type = "response")

# We need to put the threshold, in our case 0.95
my_pred = ifelse(my_pred>0.95,1,0)

# Make a contingency table to visualize the results
my_tab = table(prediction = my_pred,
               general = general_score_30$general_30)
my_tab
#             general
#prediction    0    1
#         0  397  449
#         1   35 1024

# We have:
# 397: CORRECT cancer as cancer --> 91.89815
# 1024: CORRECT control as health --> 69.51799
# 35: INCORRECT cancer as health -->  8.101852
# 449: INCORRECT health as cancer  --> 30.48201

# Calculate sensitivity and specificity 
# Sensitivity: percentage of sick people classified as sick
sensitivity = my_tab[1,1]/table(general_score_30$general_30)[1]# 0.9189815 

# Specificity: percentage of health people classified as health
specificity = my_tab[2,2]/(my_tab[1,2]+my_tab[2,2]) # 0.6951799

# Calculate AUC
auc(my_pred,general_score_30$general_30)# Area under the curve: 0.7181

# Try to check by thresholding the data
aux = general_score %>% filter(SCORE1_AVG>0.64874) # group 0 + 1st Qu.

table(aux$general_70)# 694   4

##############################################################
# OVARIAN ANALYSIS                                           #
##############################################################

# Create the file with the beast coefficients, contain: ID, Allele and freq
ovarian_variants2[,c(3,6,9)] %>% write.table("Ovarian/ovarian_coefficients.txt",
                                            row.names = F,col.names = F,quote = F)

# Calculate the ovarian score with PLINK with a maf of 0.01
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --pheno-name ovarian_70 ",
              " --score  Ovarian/ovarian_coefficients.txt ", 
              " --maf 0.01 ", 
              " --out Ovarian/ovarian_score"))

# Read the output .sscore file 
ovarian_score = fread("Ovarian/ovarian_score.sscore")
head(ovarian_score)
table(ovarian_score$ovarian_70) # 3506  295 
sum(is.na(ovarian_score$ovarian_70)) # 2768

# Filter the NAs
ovarian_score = ovarian_score %>% filter(!is.na(ovarian_70))

# Check the results of both groups 1 = control, 2 = case 
tapply(ovarian_score$SCORE1_AVG,INDEX = ovarian_score$ovarian_70,summary)

# Plot: boxplot 
boxplot(log10(ovarian_score$SCORE1_AVG)~ovarian_score$ovarian_70)

# Kruskal test 
kruskal.test(ovarian_score$SCORE1_AVG~ovarian_score$ovarian_70)

# Plot: boxplot
png("Ovarian/ovarian_boxplot_70.png")
ggplot(ovarian_score,aes(as.factor(ovarian_70),SCORE1_AVG)) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: boxplot
png("Ovarian/ovarian_boxplotlog10_70.png")
ggplot(ovarian_score,aes(as.factor(ovarian_70),log10(SCORE1_AVG))) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: histogram
png("Ovarian/ovarian_histogram_70.png")
ggplot(ovarian_score,aes(SCORE1_AVG, fill= as.factor(ovarian_70))) + 
  geom_histogram(alpha = 0.5)+
  theme_classic()
dev.off()

# Plot: histogram 
png("Ovarian/ovarian_density_70.png")
ggplot(ovarian_score,aes(log10(SCORE1_AVG), fill= as.factor(ovarian_70))) + 
  geom_density(alpha = 0.5)+
  labs(x = "Score", y = "Density log10", fill= "Groups")+
  scale_fill_discrete(name="Group",
                      breaks=c("0", "1"),
                      labels=c("Ovarian Cancer", "Health"))+
  theme_classic()
dev.off()

# As we can see almost all the cases are nead the 0 of score and the cases have higher values 

# To make the linear model (lm) we need to change the 0 and 1 
# because by default the function understand 0 as case and 1 as control 
ovarian_score$ovarian_70 = ifelse(ovarian_score$ovarian_70==2,0,1)

table(ovarian_score$ovarian_70) #  295 3506  

# Calculate the linear model with the 70% 
lm = glm(ovarian_70~SCORE1_AVG,data=ovarian_score,family = "binomial")
summary(lm)

# Make the prediction with the model, with the 70%
my_pred = predict(lm,ovarian_score,type = "response")

# We need to put the threshold, in our case 0.95
my_pred = ifelse(my_pred>0.95,1,0)

# Make a contingency table to visualize the results
my_tab = table(prediction = my_pred,
               ovarian = ovarian_score$ovarian_70)
my_tab
#             ovarian
#prediction    0      1
#         0  241     19
#         1   54   3487

# We have:
# 241: CORRECT cancer as cancer --> 81.69492
# 3487: CORRECT control as health --> 99.45807
# 54: INCORRECT cancer as health --> 18.30508
# 19: INCORRECT health as cancer  --> 0.5419281

# Calculate sensitivity and specificity 
# Sensitivity: percentage of sick people classified as sick
sensitivity = my_tab[1,1]/table(ovarian_score$ovarian_70)[1] # 0.8169492  

# Specificity: percentage of health people classified as health
specificity = my_tab[2,2]/(my_tab[1,2]+my_tab[2,2]) # 0.9945807

# Calculate AUC
auc(my_pred,ovarian_score$ovarian_70)# Area under the curve: 0.9558

# Try to check by thresholding the data
aux = ovarian_score %>% filter(SCORE1_AVG>1.35129) # group 0 + 1st Qu.

table(aux$ovarian_70)#221   3 


########################################
# Now make the prediction with the 30% #
########################################
head(ovarian_variants)
head(ovarian_score)

# Calculate the ovarian score with PLINK with a maf of 0.01
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --pheno-name ovarian_30 ",
              " --score  Ovarian/ovarian_coefficients.txt ", 
              " --maf 0.01 ", 
              " --out Ovarian/ovarian_score_30"))

# Read the output .sscore file 
ovarian_score_30 = fread("Ovarian/ovarian_score_30.sscore")
head(ovarian_score_30)
table(ovarian_score_30$ovarian_30) # 1482  129 
sum(is.na(ovarian_score_30$ovarian_30)) # 4958

# Filter the NAs
ovarian_score_30 = ovarian_score_30 %>% filter(!is.na(ovarian_30))

# Check the results of both groups 1 = control, 2 = case 
tapply(ovarian_score_30$SCORE1_AVG,INDEX = ovarian_score_30$ovarian_30,summary)

# Plot: boxplot 
boxplot(log10(ovarian_score_30$SCORE1_AVG)~ovarian_score_30$ovarian_30)

# Kruskal test 
kruskal.test(ovarian_score_30$SCORE1_AVG~ovarian_score_30$ovarian_30)

# Plot: boxplot
png("Ovarian/ovarian_boxplot_30.png")
ggplot(ovarian_score_30,aes(as.factor(ovarian_30),SCORE1_AVG)) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: boxplot
png("Ovarian/ovarian_boxplotlog10_30.png")
ggplot(ovarian_score_30,aes(as.factor(ovarian_30),log10(SCORE1_AVG))) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: histogram
png("Ovarian/ovarian_histogram_30.png")
ggplot(ovarian_score_30,aes(SCORE1_AVG, fill= as.factor(ovarian_30))) + 
  geom_histogram(alpha = 0.5)+
  theme_classic()
dev.off()

# Plot: histogram 
png("Ovarian/ovarian_density_30.png")
ggplot(ovarian_score_30,aes(log10(SCORE1_AVG), fill= as.factor(ovarian_30))) + 
  geom_density(alpha = 0.5)+
  labs(x = "Score", y = "Density log10", fill= "Groups")+
  scale_fill_discrete(name="Group",
                      breaks=c("0", "1"),
                      labels=c("Ovarian Cancer", "Health"))+
  theme_classic()
dev.off()

# To make the linear model (lm) we need to change the 0 and 1 
# because by default the function understand 0 as case and 1 as control 
ovarian_score_30$ovarian_30 = ifelse(ovarian_score_30$ovarian_30==2,0,1)

table(ovarian_score_30$ovarian_30) #  129 1482 

# Make the prediction with the model, with the 30%
my_pred = predict(lm,ovarian_score_30,type = "response")

# We need to put the threshold, in our case 0.95
my_pred = ifelse(my_pred>0.95,1,0)

# Make a contingency table to visualize the results
my_tab = table(prediction = my_pred,
               ovarian = ovarian_score_30$ovarian_30)
my_tab
#             ovarian
#prediction    0    1
#         0  109    2
#         1   20 1480

# We have:
# 109: CORRECT cancer as cancer --> 84.49612
# 1480: CORRECT control as health --> 99.86505
# 20: INCORRECT cancer as health -->  15.50388
# 2: INCORRECT health as cancer  --> 0.1349528

# Calculate sensitivity and specificity 
# Sensitivity: percentage of sick people classified as sick
sensitivity = my_tab[1,1]/table(ovarian_score_30$ovarian_30)[1] # 0.8449612 

# Specificity: percentage of health people classified as health
specificity = my_tab[2,2]/(my_tab[1,2]+my_tab[2,2]) # 0.9986505

# Calculate AUC
auc(my_pred,ovarian_score_30$ovarian_30)# Area under the curve: 0.9843

# Try to check by thresholding the data
aux = ovarian_score %>% filter(SCORE1_AVG>1.78526) # group 0 + 1st Qu.

table(aux$ovarian_70)# 210   1 


##############################################################
# BREAST ANALYSIS                                            #
##############################################################

# Create the file with the beast coefficients, contain: ID, Allele and freq
breast_variants2[,c(3,6,9)] %>% write.table("Breast/breast_coefficients.txt",
                                   row.names = F,col.names = F,quote = F)

# Calculate the breast score with PLINK with a maf of 0.01
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --pheno-name breast_70 ",
              " --score  Breast/breast_coefficients.txt ", 
              " --maf 0.01 ", 
              " --out Breast/breast_score"))

# Read the output .sscore file 
breast_score = fread("Breast/breast_score.sscore")
head(breast_score)
table(breast_score$breast_70) # 3492  584
sum(is.na(breast_score$breast_70)) # 2493

# Filter the NAs
breast_score = breast_score %>% filter(!is.na(breast_70))

# Check the results of both groups 1 = control, 2 = case 
tapply(breast_score$SCORE1_AVG,INDEX = breast_score$breast_70,summary)

# Plot: boxplot 
boxplot(log10(breast_score$SCORE1_AVG)~breast_score$breast_70)

# Kruskal test 
kruskal.test(breast_score$SCORE1_AVG~breast_score$breast_70)

# Plot: boxplot
png("Breast/breast_boxplot_70.png")
ggplot(breast_score,aes(as.factor(breast_70),SCORE1_AVG)) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: boxplot
png("Breast/breast_boxplotlog10_70.png")
ggplot(breast_score,aes(as.factor(breast_70),log10(SCORE1_AVG))) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: histogram
png("Breast/breast_histogram_70.png")
ggplot(breast_score,aes(SCORE1_AVG, fill= as.factor(breast_70))) + 
  geom_histogram(alpha = 0.5)+
  theme_classic()
dev.off()

# Plot: histogram 
png("Breast/breast_density_70.png")
ggplot(breast_score,aes(log10(SCORE1_AVG), fill= as.factor(breast_70))) + 
  geom_density(alpha = 0.5)+
  labs(x = "Score", y = "Density log10", fill= "Groups")+
  scale_fill_discrete(name="Group",
                      breaks=c("0", "1"),
                      labels=c("Breast Cancer", "Health"))+
  theme_classic()
dev.off()

# As we can see almost all the cases are nead the 0 of score and the cases have higher values 

# To make the linear model (lm) we need to change the 0 and 1 
# because by default the function understand 0 as case and 1 as control 
breast_score$breast_70 = ifelse(breast_score$breast_70==2,0,1)

table(breast_score$breast_70) #  584 3492 

# Calculate the linear model with the 70% 
lm = glm(breast_70~SCORE1_AVG,data=breast_score,family = "binomial")
summary(lm)

# Make the prediction with the model, with the 70%
my_pred = predict(lm,breast_score,type = "response")

# We need to put the threshold, in our case 0.95
my_pred = ifelse(my_pred>0.95,1,0)

# Make a contingency table to visualize the results
my_tab = table(prediction = my_pred,
                breast = breast_score$breast_70)
my_tab
#              breast
#prediction    0    1
#         0  500   76
#         1   84 3416

# We have:
# 500: CORRECT cancer as cancer --> 85.61644
# 3416: CORRECT control as health --> 97.8236
# 84: INCORRECT cancer as health --> 14.38356
# 76: INCORRECT health as cancer  --> 2.176403

# Calculate sensitivity and specificity 
# Sensitivity: percentage of sick people classified as sick
sensitivity = my_tab[1,1]/table(breast_score$breast_70)[1] # 0.8561644 

# Specificity: percentage of health people classified as health
specificity = my_tab[2,2]/(my_tab[1,2]+my_tab[2,2]) # 0.978236

# Calculate AUC
auc(my_pred,breast_score$breast_70)# Area under the curve: 0.922

# Try to check by thresholding the data
aux = breast_score %>% filter(SCORE1_AVG>0.8207) # group 0 + 1st Qu.

table(aux$breast_70)


########################################
# Now make the prediction with the 30% #
########################################
head(breast_variants)
head(breast_score)

# Calculate the breast score with PLINK with a maf of 0.01
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --pheno-name breast_30 ",
              " --score  Breast/breast_coefficients.txt ", 
              " --maf 0.01 ", 
              " --out Breast/breast_score_30"))

# Read the output .sscore file 
breast_score_30 = fread("Breast/breast_score_30.sscore")
head(breast_score_30)
table(breast_score_30$breast_30) # 1496  247 
sum(is.na(breast_score_30$breast_30)) # 4826

# Filter the NAs
breast_score_30 = breast_score_30 %>% filter(!is.na(breast_30))

# Check the results of both groups 1 = control, 2 = case 
tapply(breast_score_30$SCORE1_AVG,INDEX = breast_score_30$breast_30,summary)

# Plot: boxplot 
boxplot(log10(breast_score_30$SCORE1_AVG)~breast_score_30$breast_30)

# Kruskal test 
kruskal.test(breast_score_30$SCORE1_AVG~breast_score_30$breast_30)

# Plot: boxplot
png("Breast/breast_boxplot_30.png")
ggplot(breast_score_30,aes(as.factor(breast_30),SCORE1_AVG)) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: boxplot
png("Breast/breast_boxplotlog10_30.png")
ggplot(breast_score_30,aes(as.factor(breast_30),log10(SCORE1_AVG))) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: histogram
png("Breast/breast_histogram_30.png")
ggplot(breast_score_30,aes(SCORE1_AVG, fill= as.factor(breast_30))) + 
  geom_histogram(alpha = 0.5)+
  theme_classic()
dev.off()

# Plot: histogram 
png("Breast/breast_density_30.png")
ggplot(breast_score_30,aes(log10(SCORE1_AVG), fill= as.factor(breast_30))) + 
  geom_density(alpha = 0.5)+
  labs(x = "Score", y = "Density log10", fill= "Groups")+
  scale_fill_discrete(name="Group",
                      breaks=c("0", "1"),
                      labels=c("Breast Cancer", "Health"))+
  theme_classic()
dev.off()

# To make the linear model (lm) we need to change the 0 and 1 
# because by default the function understand 0 as case and 1 as control 
breast_score_30$breast_30 = ifelse(breast_score_30$breast_30==2,0,1)

table(breast_score_30$breast_30) #  247 1496

# Make the prediction with the model, with the 30%
my_pred = predict(lm,breast_score_30,type = "response")

# We need to put the threshold, in our case 0.95
my_pred = ifelse(my_pred>0.95,1,0)

# Make a contingency table to visualize the results
my_tab = table(prediction = my_pred,
               breast = breast_score_30$breast_30)
my_tab
#             breast
#prediction    0    1
#         0  217   30
#         1   30 1466

# We have:
# 217: CORRECT cancer as cancer --> 87.85425
# 1466: CORRECT control as health --> 97.99465
# 30: INCORRECT cancer as health -->  12.14575
# 30: INCORRECT health as cancer  --> 2.005348

# Calculate sensitivity and specificity 
# Sensitivity: percentage of sick people classified as sick
sensitivity = my_tab[1,1]/table(breast_score_30$breast_30)[1] # 0.8785425 

# Specificity: percentage of health people classified as health
specificity = my_tab[2,2]/(my_tab[1,2]+my_tab[2,2]) # 0.9799465

# Calculate AUC
auc(my_pred,breast_score_30$breast_30)# Area under the curve: 0.9292

# Try to check by thresholding the data
aux = breast_score %>% filter(SCORE1_AVG>0.95195) # group 0 + 1st Qu.

table(aux$breast_70)# 414   3 


##############################################################
# COLON ANALYSIS                                            #
##############################################################

# Create the file with the beast coefficients, contain: ID, Allele and freq
colon_variants2[,c(3,6,9)] %>% write.table("Colon/colon_coefficients.txt",
                                            row.names = F,col.names = F,quote = F)

# Calculate the colon score with PLINK with a maf of 0.01
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --pheno-name colon_70 ",
              " --score  Colon/colon_coefficients.txt ", 
              " --maf 0.01 ", 
              " --out Colon/colon_score"))

# Read the output .sscore file 
colon_score = fread("Colon/colon_score.sscore")
head(colon_score)
table(colon_score$colon_70) # 3489  219
sum(is.na(colon_score$colon_70)) # 2861

# Filter the NAs
colon_score = colon_score %>% filter(!is.na(colon_70))

# Check the results of both groups 1 = control, 2 = case 
tapply(colon_score$SCORE1_AVG,INDEX = colon_score$colon_70,summary)

# Plot: boxplot 
boxplot(log10(colon_score$SCORE1_AVG)~colon_score$colon_70)

# Kruskal test 
kruskal.test(colon_score$SCORE1_AVG~colon_score$colon_70)

# Plot: boxplot
png("Colon/colon_boxplot_70.png")
ggplot(colon_score,aes(as.factor(colon_70),SCORE1_AVG)) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: boxplot
png("Colon/colon_boxplotlog10_70.png")
ggplot(colon_score,aes(as.factor(colon_70),log10(SCORE1_AVG))) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: histogram
png("Colon/colon_histogram_70.png")
ggplot(colon_score,aes(SCORE1_AVG, fill= as.factor(colon_70))) + 
  geom_histogram(alpha = 0.5)+
  theme_classic()
dev.off()

# Plot: histogram 
png("Colon/colon_density_70.png")
ggplot(colon_score,aes(log10(SCORE1_AVG), fill= as.factor(colon_70))) + 
  geom_density(alpha = 0.5)+
  labs(x = "Score", y = "Density log10", fill= "Groups")+
  scale_fill_discrete(name="Group",
                      breaks=c("0", "1"),
                      labels=c("Colon Cancer", "Health"))+
  theme_classic()
dev.off()

# As we can see almost all the cases are nead the 0 of score and the cases have higher values 

# To make the linear model (lm) we need to change the 0 and 1 
# because by default the function understand 0 as case and 1 as control 
colon_score$colon_70 = ifelse(colon_score$colon_70==2,0,1)

table(colon_score$colon_70) #  219 3489 

# Calculate the linear model with the 70% 
lm = glm(colon_70~SCORE1_AVG,data=colon_score,family = "binomial")
summary(lm)

# Make the prediction with the model, with the 70%
my_pred = predict(lm,colon_score,type = "response")

# We need to put the threshold, in our case 0.95
my_pred = ifelse(my_pred>0.95,1,0)

# Make a contingency table to visualize the results
my_tab = table(prediction = my_pred,
               colon = colon_score$colon_70)
my_tab
#           colon
#prediction    0    1
#         0  179   11
#         1   40 3478

# We have:
# 179: CORRECT cancer as cancer --> 81.73516
# 3478: CORRECT control as health --> 99.68472
# 40: INCORRECT cancer as health --> 18.26484
# 11: INCORRECT health as cancer  --> 0.3152766

# Calculate sensitivity and specificity 
# Sensitivity: percentage of sick people classified as sick
sensitivity = my_tab[1,1]/table(colon_score$colon_70)[1] # 0.8173516 

# Specificity: percentage of health people classified as health
specificity = my_tab[2,2]/(my_tab[1,2]+my_tab[2,2]) # 0.9968472

# Calculate AUC
auc(my_pred,colon_score$colon_70)# Area under the curve: 0.9654

# Try to check by thresholding the data
aux = colon_score %>% filter(SCORE1_AVG>4.14404) # group 0 + 1st Qu.

table(aux$colon_70) # 164 


########################################
# Now make the prediction with the 30% #
########################################
head(colon_variants)
head(colon_score)

# Calculate the colon score with PLINK with a maf of 0.01
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --pheno-name colon_30 ",
              " --score  Colon/colon_coefficients.txt ", 
              " --maf 0.01 ", 
              " --out Colon/colon_score_30"))

# Read the output .sscore file 
colon_score_30 = fread("Colon/colon_score_30.sscore")
head(colon_score_30)
table(colon_score_30$colon_30) # 1496  247 
sum(is.na(colon_score_30$colon_30)) # 4826

# Filter the NAs
colon_score_30 = colon_score_30 %>% filter(!is.na(colon_30))

# Check the results of both groups 1 = control, 2 = case 
tapply(colon_score_30$SCORE1_AVG,INDEX = colon_score_30$colon_30,summary)

# Plot: boxplot 
boxplot(log10(colon_score_30$SCORE1_AVG)~colon_score_30$colon_30)

# Kruskal test 
kruskal.test(colon_score_30$SCORE1_AVG~colon_score_30$colon_30)

# Plot: boxplot
png("Colon/colon_boxplot_30.png")
ggplot(colon_score_30,aes(as.factor(colon_30),SCORE1_AVG)) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: boxplot
png("Colon/colon_boxplotlog10_30.png")
ggplot(colon_score_30,aes(as.factor(colon_30),log10(SCORE1_AVG))) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8)+
  labs(x = "Groups", y = "Score", fill= "Groups", title = "Score boxplot by groups: 2 = Case and 1 = Control")+
  theme_classic()
dev.off()

# Plot: histogram
png("Colon/colon_histogram_30.png")
ggplot(colon_score_30,aes(SCORE1_AVG, fill= as.factor(colon_30))) + 
  geom_histogram(alpha = 0.5)+
  theme_classic()
dev.off()

# Plot: histogram 
png("Colon/colon_density_30.png")
ggplot(colon_score_30,aes(log10(SCORE1_AVG), fill= as.factor(colon_30))) + 
  geom_density(alpha = 0.5)+
  labs(x = "Score", y = "Density log10", fill= "Groups")+
  scale_fill_discrete(name="Group",
                      breaks=c("0", "1"),
                      labels=c("Colon Cancer", "Health"))+
  theme_classic()
dev.off()

# To make the linear model (lm) we need to change the 0 and 1 
# because by default the function understand 0 as case and 1 as control 
colon_score_30$colon_30 = ifelse(colon_score_30$colon_30==2,0,1)

table(colon_score_30$colon_30) #  100 1499 

# Make the prediction with the model, with the 30%
my_pred = predict(lm,colon_score_30,type = "response")

# We need to put the threshold, in our case 0.95
my_pred = ifelse(my_pred>0.95,1,0)

# Make a contingency table to visualize the results
my_tab = table(prediction = my_pred,
               colon = colon_score_30$colon_30)
my_tab
#               colon
#prediction    0    1
#         0   71    7
#         1   29 1492

# We have:
# 71: CORRECT cancer as cancer --> 71
# 1492: CORRECT control as health --> 99.53302
# 29: INCORRECT cancer as health -->  29
# 7: INCORRECT health as cancer  --> 0.466978

# Calculate sensitivity and specificity 
# Sensitivity: percentage of sick people classified as sick
sensitivity = my_tab[1,1]/table(colon_score_30$colon_30)[1] #0.71

# Specificity: percentage of health people classified as health
specificity = my_tab[2,2]/(my_tab[1,2]+my_tab[2,2]) # 0.9953302

# Calculate AUC
auc(my_pred,colon_score_30$colon_30)# Area under the curve: 0.9456

# Try to check by thresholding the data
aux = colon_score %>% filter(SCORE1_AVG>0.41060) # group 0 + 1st Qu.

table(aux$colon_70)# 184  23 
