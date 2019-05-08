######################################################################################
# 08 PRS advanced                                                                    #
######################################################################################

# Set the working directory
setwd("/imppc/labs/dnalab/share/PRS/08_PRS_advanced")

# Libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(pROC)
library(caret)
library(e1071)
library(xlsx)
library(grid)


# Create the table with all the posible combinations
prs_statistics = expand.grid(MAF_cases = c(0.001,0.01,0.05),
                             p_value_threshold = c(0.05,0.01,1e-3,1e-4,1e-5,1e-6,5e-8),
                             LD_pruning = c("YES","NO"),
                             classification_threshold = c(0.05,0.07,0.10, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55))


head(prs_statistics)

# Add the columns we will have to fill
prs_statistics$N_snps = 0
prs_statistics$sensitivity = 0
prs_statistics$specificity = 0
prs_statistics$AUC = 0


# generate list of variants in LD

#dir.create("outputs")

system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --indep-pairwise 50 5 0.2  ",
              " --out /imppc/labs/dnalab/share/PRS/08_PRS_advanced/outputs/LD"))

ld_in = fread("outputs/LD.prune.in",header = F)

# Check the number
dim(ld_in) # 8689


## LOOP
for(i in 1:nrow(prs_statistics)){ # For each row of the table (each row is a variables combination)
  # Take the variables values for this combination
  maf = prs_statistics$MAF_cases[i]
  p_value = prs_statistics$p_value_threshold[i]
  ld = prs_statistics$LD_pruning[i]
  threshold = prs_statistics$classification_threshold[i]
  
  # Perform the GWAS
  system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
                " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
                " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
                " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
                " --pheno-name pheno_70 ",
                " --glm firth cols=chrom,pos,ref,alt,firth,test,nobs,orbeta,se,ci,tz,p,a1countcc,a1freqcc,totallelecc ",  
                " --out /imppc/labs/dnalab/share/PRS/08_PRS_advanced/outputs/gwas"))
  # Read the GWAS
  gwas = fread("/imppc/labs/dnalab/share/PRS/08_PRS_advanced/outputs/gwas.pheno_70.glm.firth")
  
  # Apply the filters 
  gwas = gwas %>% filter(A1_CASE_FREQ>maf & P<p_value)
  
  # Condition of the LD 
  if(ld=="YES"){
    gwas = gwas[gwas$ID %in% ld_in$V1,] # If the analysis is with LD we filter the variants in the ld_in 
  }
  
  # Select the GWAS data we need to calculate the score
  gwas[,c(3,6,15)]  %>% write.table(paste0("outputs/snps_score_",maf,"_",p_value,"_",ld),
                              row.names=F,quote=F)
  
  # Calculate the score for 70%
  system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
                " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
                " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
                " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
                " --pheno-name pheno_70 ",
                " --score  outputs/snps_score_",maf,"_",p_value,"_",ld, 
                " --out outputs/prs_",maf,"_",p_value,"_",ld,"_",threshold,"_70"))
  
  # Read the score and change column name
  general_score = fread(paste0("outputs/prs_",maf,"_",p_value,"_",ld,"_",threshold,"_70.sscore"))
  colnames(general_score)[3] = "general_70"
  
  # Delete na values
  general_score = general_score %>% filter(!is.na(general_70))
  
  # Change the 1 and 2 values for 0 and 1 (the ones we need to develop the model)
  general_score$general_70 = ifelse(general_score$general_70==1,0,1)
  
  # indicate the characteristics of the train function
  ctrl <- trainControl(method = "repeatedcv", number = 5,
                       savePredictions = TRUE)
  
  # Calculate the linear model with the 70% 
  lm_general = train(general_70~SCORE1_AVG,data=general_score,method="glm", family="binomial",
                     trControl = ctrl)
  
  # Calculate the score for 30%
  system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
                " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
                " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
                " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
                " --pheno-name pheno_30 ",
                " --score  outputs/snps_score_",maf,"_",p_value,"_",ld, 
                " --out outputs/prs_",maf,"_",p_value,"_",ld,"_",threshold,"_30"))
  
  # Read the scores file and change the column name 
  general_score_30 = fread(paste0("outputs/prs_",maf,"_",p_value,"_",ld,"_",threshold,"_30.sscore"))
  colnames(general_score_30)[3] = "general_30"
  
  # Filter the NAs
  general_score_30 = general_score_30 %>% filter(!is.na(general_30))
  
    # To make the linear model (lm) we need to change the 0 and 1 
  # because by default the function understand 0 as case and 1 as control 
  general_score_30$general_30 = ifelse(general_score_30$general_30==2,1,0)
  
  # Make the prediction with the model, with the 30%
  my_pred_general = predict(lm_general,general_score_30,type = "raw")
  
  # We need to put the threshold, in our case we use different thresholds
  my_pred_general = ifelse(my_pred_general>threshold,1,0)
  
  # Make a contingency table to visualize the results
  my_tab = table(prediction = my_pred_general,
                 general = general_score_30$general_30)
  
  # Sensitivity: percentage of sick people classified as sick
  sensitivity = my_tab[2,2]/table(general_score_30$general_30)[2]
  
  # Specificity: percentage of health people classified as health
  specificity = my_tab[1,1]/table(general_score_30$general_30)[1]
  
  # Calculate AUC
  my_auc = auc(my_pred_general,general_score_30$general_30)
  
  # Save all this values we calculate in the table
  prs_statistics$N_snps[i] = nrow(gwas)
  prs_statistics$sensitivity[i] = sensitivity
  prs_statistics$specificity[i] = specificity
  prs_statistics$AUC[i] = my_auc
  
}

# Check the table
prs_statistics

# Look only at the LD_pruning analysis combinations
prs_statistics %>% filter(LD_pruning == "YES") %>% arrange(desc(sensitivity))

# For maf 0.05
best005 = prs_statistics %>% 
              filter(MAF_cases == 0.05 & LD_pruning == "YES" & specificity > 0.80) %>% 
              arrange(desc(sensitivity)) %>% 
              head()

# For maf 0.01
best001 = prs_statistics %>% 
              filter(MAF_cases == 0.01 & LD_pruning == "YES" & specificity > 0.80) %>% 
              arrange(desc(sensitivity)) %>% 
              head()

# For maf 0.001
best0001 = prs_statistics %>% 
              filter(MAF_cases == 0.001 & LD_pruning == "YES" & specificity > 0.80) %>% 
              arrange(desc(sensitivity)) %>% 
              head()

best_prs = rbind(best005, best001, best0001)

# Save the information:
write.xlsx(x = best_prs, file = "best_prs.xlsx",
           sheetName = "PRS_statistics", row.names = FALSE)

best_prs %>% write.table("best_prs.txt",
                                  row.names = F,col.names = F,quote = F)

prs_statistics %>% write.table("prs_statistics.txt",
                               row.names = F,col.names = F,quote = F)

#############################################################################################
# BEST PRS MODEL FOR MAF 0.05                                                               #
#############################################################################################
# p_value_threshold = 1e-06, Pruning YES, Classification_threshold = 0.07

#Read 70 scores
general_score_70 = fread(paste0("outputs/prs_0.05_1e-06_YES_0.07_70.sscore"))
colnames(general_score_70)[3] = "general"
general_score_70 = general_score_70 %>% filter(!is.na(general))

#Read 30 scores
general_score_30 = fread(paste0("outputs/prs_0.05_1e-06_YES_0.07_30.sscore"))
colnames(general_score_30)[3] = "general"
general_score_30 = general_score_30 %>% filter(!is.na(general))

head(general_score_30)
head(general_score_70)

# Merge the both scores
general_score = rbind(general_score_70,general_score_30)
head(general_score)
table(general_score$general)

# Order the scores by SCORE1_AVG
general_score = general_score %>% arrange(SCORE1_AVG)

# Create the prs_percentile variable
general_score$prs_percentile = 99

# Create the percentiles first number
my_seq = seq(1,nrow(general_score),round(nrow(general_score)*0.01))

# Fill the percentiles variable 
cont = 1
for(i in 1:(length(my_seq)-1)){
  
  general_score$prs_percentile[my_seq[i]:(my_seq[i+1]-1)] = cont
  
  cont = cont+1  
}

# Check the percentiles N
table(general_score$prs_percentile)

# Create the prs_percentiles dataframe 
prs_percentile = data.frame(percentile=seq(1,99,1))

# Create the prevalence variable
prs_percentile$prevalence = 0
prs_percentile$meanscore = 0

# Fill the prevalence
for(i in 1:nrow(prs_percentile)){
  # Obtain the percentile information
  aux = general_score %>% filter(prs_percentile==i)
  # Case and control information
  aux2 = prop.table(table(aux$general))
  # Obtain the average score of the percentile
  prs_percentile$meanscore[i] = mean(aux$SCORE1_AVG)
  # Fill the prevalence
  if (dim(aux2) == 1){
    if (aux$general[1]== 1){
      prs_percentile$prevalence[i] = 0
    }
    if (aux$general[1]== 2){
      prs_percentile$prevalence[i] = 100
    }
  }
  if (dim(aux2) == 2){
    prs_percentile$prevalence[i] = aux2[2]*100 
  }
   
}

# Add information for colour
prs_percentile$colorthreshold = 0
prs_percentile[81:99,]$colorthreshold = 1


# Check 
prs_percentile

# Boxplot
boxplot(log(general_score$SCORE1_AVG)~general_score$general)

# Plot the prs_percentile 
png("graphs/005/005_prevalence_percentile.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
p <-  ggplot(prs_percentile, aes(percentile, prevalence)) +
      geom_point(color = "dodgerblue3") +
      labs(title = "Prevalence of cancer by percentile PRS maf 0.05", 
           subtitle = "All samples, GCAT and Pandora")+
      theme_classic() 
p + theme(
    plot.title = element_text(color = "dodgerblue4", size = 153 , face = "bold"),
    plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
    axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
    axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")
)
dev.off()

# Plot the prs_percentile  with colour by score
png("graphs/005/005_prevalence_percentile_colorbyscore.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
  q <- ggplot(prs_percentile, aes(percentile, prevalence, color = (meanscore))) +
        geom_point() +
        labs(title = "Prevalence of cancer by percentile PRS maf 0.05", 
             subtitle = "All samples, GCAT and Pandora \nColour by increase in score") +
        scale_fill_discrete(name="Mean Score") +
        theme_classic()
  q + theme(
    plot.title = element_text(color = "dodgerblue4", size = 153 , face = "bold"),
    plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
    axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
    axis.title.y = element_text(color="dodgerblue4", size=10, face="italic"),
    legend.title=element_blank())
dev.off()

# Plot the prs_percentile  with colour by score WITH LOG10
png("graphs/005/005_prevalence_percentile_colorbyscore_log10.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
a <- ggplot(prs_percentile, aes(percentile, prevalence, color = log10(meanscore))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.05", 
       subtitle = "All samples, GCAT and Pandora \nColour by increase in score") +
  scale_fill_discrete(name="Mean Score log10") +
  theme_classic()
a + theme(
  plot.title = element_text(color = "dodgerblue4", size = 153 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic"),
  legend.title=element_blank())
dev.off()

# Plot the prs_percentile  with colour by score WITH LOG10
png("graphs/005/005_prevalence_percentile_colorbythreshold.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
b <- ggplot(prs_percentile, aes(x = percentile, y = prevalence, color = as.factor(colorthreshold))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.05", 
       subtitle = "All samples, GCAT and Pandora \nColour by higher percentage of cases or control") +
  theme_classic()
b + scale_color_manual(values = c("lightskyblue", "dodgerblue4"),
                      name  ="",
                      breaks=c("0", "1"),
                      labels=c("High control \npercentage \n", "High case \npercentage")) +
  theme(
    plot.title = element_text(color = "dodgerblue4", size = 153 , face = "bold"),
    plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
    axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
    axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_hline(yintercept=50, linetype="dashed", color = "black") +
  geom_vline(xintercept=83.5, linetype="dashed", color = "black")
dev.off()



# Now we want to visualize GCAT and Pandora separately
head(general_score)
table(general_score$general)

general_score_gcat = general_score[which(general_score$general == 1),] 
general_score_pandora = general_score[which(general_score$general == 2),]


head(general_score_gcat)
dim(general_score_gcat)


# Create the prs_percentile variable
general_score_gcat$prs_percentile = 99
general_score_pandora$prs_percentile = 99

# Create the percentiles first number
my_seq_gcat = seq(1,nrow(general_score_gcat),round(nrow(general_score_gcat)*0.01))
my_seq_pandora = seq(1,nrow(general_score_pandora),round(nrow(general_score_pandora)*0.01))

# Fill the percentiles variable 
cont = 1
for(i in 1:(length(my_seq_gcat)-1)){
  
  general_score_gcat$prs_percentile[my_seq_gcat[i]:(my_seq_gcat[i+1]-1)] = cont
  
  cont = cont+1  
}
cont = 1
for(i in 1:(length(my_seq_pandora)-1)){
  
  general_score_pandora$prs_percentile[my_seq_pandora[i]:(my_seq_pandora[i+1]-1)] = cont
  
  cont = cont+1  
}


# Check the percentiles N
table(general_score_gcat$prs_percentile)
table(general_score_pandora$prs_percentile)
head(general_score_gcat)
head(general_score_pandora)

# Create the prs_percentiles dataframe 
prs_percentile_gcat = data.frame(percentile=seq(1,99,1))
prs_percentile_pandora = data.frame(percentile=seq(1,99,1))

# Create the prevalence variable
prs_percentile_gcat$meanscore = 0
prs_percentile_pandora$meanscore = 0

# Fill the prevalence GCAT and Pandora
for(i in 1:nrow(prs_percentile_gcat)){
  # Obtain the percentile information
  aux_G = general_score_gcat %>% filter(prs_percentile==i)
  aux_P = general_score_pandora %>% filter(prs_percentile==i)
  # Obtain the average score of the percentile
  prs_percentile_gcat$meanscore[i] = mean(aux_G$SCORE1_AVG)
  prs_percentile_pandora$meanscore[i] = mean(aux_P$SCORE1_AVG)
}

# Check 
prs_percentile_gcat
prs_percentile_pandora

# Change column names and merge the both datasets
colnames(prs_percentile_gcat)[2] = "meanscore_GCAT"
colnames(prs_percentile_pandora)[2] = "meanscore_Pandora"
prs_percentile_gcat_pandora = merge(prs_percentile_gcat, prs_percentile_pandora)
head(prs_percentile_gcat_pandora)


png("graphs/005/005_scoreincrease_bygroup.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
grob <- grobTree(textGrob("Score = 1 \npercentile = 19", x=0.05,  y=0.15, hjust=0,
                          gp=gpar(col="red", fontsize=8, fontface="italic")))
grob2 <- grobTree(textGrob("948 GCAT \n261 Pandora \n19%", x=0.05,  y=0.8, hjust=0,
                          gp=gpar(col="black", fontsize=8, fontface="italic")))
grob3 <- grobTree(textGrob("4040 GCAT \n1114 Pandora \n81%", x=0.25,  y=0.8, hjust=0,
                          gp=gpar(col="black", fontsize=8, fontface="italic")))
e <- ggplot(prs_percentile_gcat_pandora) +
  geom_line(aes(y= prs_percentile_gcat_pandora$meanscore_GCAT,
                x = prs_percentile_gcat_pandora$percentile),
            colour= "lightskyblue") +
  geom_line(aes(y= prs_percentile_gcat_pandora$meanscore_Pandora,
                x = prs_percentile_gcat_pandora$percentile),
            colour = "dodgerblue4") +
  labs(title = "Increase in Score by group PRS maf 0.05", 
       subtitle = "GCAT as ligthblue and Pandora as darkblue")+
  xlab("percentile")+
  ylab("PRS") +
  theme_classic()
e + theme(
    plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
    plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
    axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
    axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  geom_vline(xintercept=19, linetype="dashed", color = "red") +
  annotation_custom(grob) + annotation_custom(grob2) + annotation_custom(grob3)
dev.off()

## Histogram plot

ggplot(general_score,aes(log10(SCORE1_AVG), fill = as.factor(general)))+
  geom_density(alpha = 0.5)

ggplot(general_score,aes(log10(SCORE1_AVG), colour = as.factor(general)))+
  geom_freqpoly(binwidth = 0.05)


#############################################################################################
# BEST PRS MODEL FOR MAF 0.01                                                               #
#############################################################################################
# p_value_threshold = 1e-04, Pruning YES, Classification_threshold = 0.05

#Read 70 scores
general_score_70 = fread(paste0("outputs/prs_0.01_1e-04_YES_0.05_70.sscore"))
colnames(general_score_70)[3] = "general"
general_score_70 = general_score_70 %>% filter(!is.na(general))

#Read 30 scores
general_score_30 = fread(paste0("outputs/prs_0.01_1e-04_YES_0.05_30.sscore"))
colnames(general_score_30)[3] = "general"
general_score_30 = general_score_30 %>% filter(!is.na(general))

head(general_score_30)
head(general_score_70)

# Merge the both scores
general_score = rbind(general_score_70,general_score_30)
head(general_score)
table(general_score$general)

# Order the scores by SCORE1_AVG
general_score = general_score %>% arrange(SCORE1_AVG)

# Create the prs_percentile variable
general_score$prs_percentile = 99

# Create the percentiles first number
my_seq = seq(1,nrow(general_score),round(nrow(general_score)*0.01))

# Fill the percentiles variable 
cont = 1
for(i in 1:(length(my_seq)-1)){
  
  general_score$prs_percentile[my_seq[i]:(my_seq[i+1]-1)] = cont
  
  cont = cont+1  
}

# Check the percentiles N
table(general_score$prs_percentile)

# Create the prs_percentiles dataframe 
prs_percentile = data.frame(percentile=seq(1,99,1))

# Create the prevalence variable
prs_percentile$prevalence = 0
prs_percentile$meanscore = 0

# Fill the prevalence
for(i in 1:nrow(prs_percentile)){
  # Obtain the percentile information
  aux = general_score %>% filter(prs_percentile==i)
  # Case and control information
  aux2 = prop.table(table(aux$general))
  # Obtain the average score of the percentile
  prs_percentile$meanscore[i] = mean(aux$SCORE1_AVG)
  # Fill the prevalence
  if (dim(aux2) == 1){
    if (aux$general[1]== 1){
      prs_percentile$prevalence[i] = 0
    }
    if (aux$general[1]== 2){
      prs_percentile$prevalence[i] = 100
    }
  }
  if (dim(aux2) == 2){
    prs_percentile$prevalence[i] = aux2[2]*100 
  }
  
}

# Add information for colour
prs_percentile$colorthreshold = 0
prs_percentile[81:99,]$colorthreshold = 1


# Check 
prs_percentile

# Boxplot
boxplot(log(general_score$SCORE1_AVG)~general_score$general)

# Plot the prs_percentile 
png("graphs/001/001_prevalence_percentile.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
p <-  ggplot(prs_percentile, aes(percentile, prevalence)) +
  geom_point(color = "dodgerblue3") +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.01", 
       subtitle = "All samples, GCAT and Pandora")+
  theme_classic() 
p + theme(
  plot.title = element_text(color = "dodgerblue4", size = 153 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")
)
dev.off()

# Plot the prs_percentile  with colour by score
png("graphs/001/001_prevalence_percentile_colorbyscore.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
q <- ggplot(prs_percentile, aes(percentile, prevalence, color = (meanscore))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.01", 
       subtitle = "All samples, GCAT and Pandora \nColour by increase in score") +
  scale_fill_discrete(name="Mean Score") +
  theme_classic()
q + theme(
  plot.title = element_text(color = "dodgerblue4", size = 153 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic"),
  legend.title=element_blank())
dev.off()

# Plot the prs_percentile  with colour by score WITH LOG10
png("graphs/001/001_prevalence_percentile_colorbyscore_log10.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
a <- ggplot(prs_percentile, aes(percentile, prevalence, color = log10(meanscore))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.01", 
       subtitle = "All samples, GCAT and Pandora \nColour by increase in score") +
  scale_fill_discrete(name="Mean Score log10") +
  theme_classic()
a + theme(
  plot.title = element_text(color = "dodgerblue4", size = 153 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic"),
  legend.title=element_blank())
dev.off()

# Plot the prs_percentile  with colour by score WITH LOG10
png("graphs/001/001_prevalence_percentile_colorbythreshold.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
b <- ggplot(prs_percentile, aes(x = percentile, y = prevalence, color = as.factor(colorthreshold))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.01", 
       subtitle = "All samples, GCAT and Pandora \nColour by higher percentage of cases or control") +
  theme_classic()
b + scale_color_manual(values = c("lightskyblue", "dodgerblue4"),
                       name  ="",
                       breaks=c("0", "1"),
                       labels=c("High control \npercentage \n", "High case \npercentage")) +
  theme(
    plot.title = element_text(color = "dodgerblue4", size = 153 , face = "bold"),
    plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
    axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
    axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_hline(yintercept=50, linetype="dashed", color = "black") +
  geom_vline(xintercept=83.5, linetype="dashed", color = "black")
dev.off()



# Now we want to visualize GCAT and Pandora separately
head(general_score)
table(general_score$general)

general_score_gcat = general_score[which(general_score$general == 1),] 
general_score_pandora = general_score[which(general_score$general == 2),]


head(general_score_gcat)
dim(general_score_gcat)


# Create the prs_percentile variable
general_score_gcat$prs_percentile = 99
general_score_pandora$prs_percentile = 99

# Create the percentiles first number
my_seq_gcat = seq(1,nrow(general_score_gcat),round(nrow(general_score_gcat)*0.01))
my_seq_pandora = seq(1,nrow(general_score_pandora),round(nrow(general_score_pandora)*0.01))

# Fill the percentiles variable 
cont = 1
for(i in 1:(length(my_seq_gcat)-1)){
  
  general_score_gcat$prs_percentile[my_seq_gcat[i]:(my_seq_gcat[i+1]-1)] = cont
  
  cont = cont+1  
}
cont = 1
for(i in 1:(length(my_seq_pandora)-1)){
  
  general_score_pandora$prs_percentile[my_seq_pandora[i]:(my_seq_pandora[i+1]-1)] = cont
  
  cont = cont+1  
}


# Check the percentiles N
table(general_score_gcat$prs_percentile)
table(general_score_pandora$prs_percentile)
head(general_score_gcat)
head(general_score_pandora)

# Create the prs_percentiles dataframe 
prs_percentile_gcat = data.frame(percentile=seq(1,99,1))
prs_percentile_pandora = data.frame(percentile=seq(1,99,1))

# Create the prevalence variable
prs_percentile_gcat$meanscore = 0
prs_percentile_pandora$meanscore = 0

# Fill the prevalence GCAT and Pandora
for(i in 1:nrow(prs_percentile_gcat)){
  # Obtain the percentile information
  aux_G = general_score_gcat %>% filter(prs_percentile==i)
  aux_P = general_score_pandora %>% filter(prs_percentile==i)
  # Obtain the average score of the percentile
  prs_percentile_gcat$meanscore[i] = mean(aux_G$SCORE1_AVG)
  prs_percentile_pandora$meanscore[i] = mean(aux_P$SCORE1_AVG)
}

# Check 
prs_percentile_gcat
prs_percentile_pandora

# Change column names and merge the both datasets
colnames(prs_percentile_gcat)[2] = "meanscore_GCAT"
colnames(prs_percentile_pandora)[2] = "meanscore_Pandora"
prs_percentile_gcat_pandora = merge(prs_percentile_gcat, prs_percentile_pandora)
head(prs_percentile_gcat_pandora)


png("graphs/001/001_scoreincrease_bygroup.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
grob <- grobTree(textGrob("Score = 0.6 \npercentile = 15", x=0.005,  y=0.15, hjust=0,
                          gp=gpar(col="red", fontsize=8, fontface="italic")))
grob2 <- grobTree(textGrob("748 GCAT \n206 Pandora \n15%", x=0.02,  y=0.8, hjust=0,
                           gp=gpar(col="black", fontsize=8, fontface="italic")))
grob3 <- grobTree(textGrob("4240 GCAT \n1169 Pandora \n85%", x=0.2,  y=0.8, hjust=0,
                           gp=gpar(col="black", fontsize=8, fontface="italic")))
e <- ggplot(prs_percentile_gcat_pandora) +
  geom_line(aes(y= prs_percentile_gcat_pandora$meanscore_GCAT,
                x = prs_percentile_gcat_pandora$percentile),
            colour= "lightskyblue") +
  geom_line(aes(y= prs_percentile_gcat_pandora$meanscore_Pandora,
                x = prs_percentile_gcat_pandora$percentile),
            colour = "dodgerblue4") +
  labs(title = "Increase in Score by group PRS maf 0.01", 
       subtitle = "GCAT as ligthblue and Pandora as darkblue")+
  xlab("percentile")+
  ylab("PRS") +
  theme_classic()
e + theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_hline(yintercept=0.6, linetype="dashed", color = "red") +
  geom_vline(xintercept=15, linetype="dashed", color = "red") +
  annotation_custom(grob) + annotation_custom(grob2) + annotation_custom(grob3)
dev.off()



#############################################################################################
# BEST PRS MODEL FOR MAF 0.001                                                               #
#############################################################################################
# p_value_threshold = 0.05, Pruning YES, Classification_threshold = 0.05


#Read 70 scores
general_score_70 = fread(paste0("outputs/prs_0.001_0.05_YES_0.05_70.sscore"))
colnames(general_score_70)[3] = "general"
general_score_70 = general_score_70 %>% filter(!is.na(general))

#Read 30 scores
general_score_30 = fread(paste0("outputs/prs_0.001_0.05_YES_0.05_30.sscore"))
colnames(general_score_30)[3] = "general"
general_score_30 = general_score_30 %>% filter(!is.na(general))

head(general_score_30)
head(general_score_70)

# Merge the both scores
general_score = rbind(general_score_70,general_score_30)
head(general_score)
table(general_score$general)

# Order the scores by SCORE1_AVG
general_score = general_score %>% arrange(SCORE1_AVG)

# Create the prs_percentile variable
general_score$prs_percentile = 99

# Create the percentiles first number
my_seq = seq(1,nrow(general_score),round(nrow(general_score)*0.01))

# Fill the percentiles variable 
cont = 1
for(i in 1:(length(my_seq)-1)){
  
  general_score$prs_percentile[my_seq[i]:(my_seq[i+1]-1)] = cont
  
  cont = cont+1  
}

# Check the percentiles N
table(general_score$prs_percentile)

# Create the prs_percentiles dataframe 
prs_percentile = data.frame(percentile=seq(1,99,1))

# Create the prevalence variable
prs_percentile$prevalence = 0
prs_percentile$meanscore = 0

# Fill the prevalence
for(i in 1:nrow(prs_percentile)){
  # Obtain the percentile information
  aux = general_score %>% filter(prs_percentile==i)
  # Case and control information
  aux2 = prop.table(table(aux$general))
  # Obtain the average score of the percentile
  prs_percentile$meanscore[i] = mean(aux$SCORE1_AVG)
  # Fill the prevalence
  if (dim(aux2) == 1){
    if (aux$general[1]== 1){
      prs_percentile$prevalence[i] = 0
    }
    if (aux$general[1]== 2){
      prs_percentile$prevalence[i] = 100
    }
  }
  if (dim(aux2) == 2){
    prs_percentile$prevalence[i] = aux2[2]*100 
  }
  
}

# Add information for colour
prs_percentile$colorthreshold = 0
prs_percentile[81:99,]$colorthreshold = 1


# Check 
prs_percentile

# Boxplot
boxplot(log(general_score$SCORE1_AVG)~general_score$general)

# Plot the prs_percentile 
png("graphs/0001/0001_prevalence_percentile.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
p <-  ggplot(prs_percentile, aes(percentile, prevalence)) +
  geom_point(color = "dodgerblue3") +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.001", 
       subtitle = "All samples, GCAT and Pandora")+
  theme_classic() 
p + theme(
  plot.title = element_text(color = "dodgerblue4", size = 153 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")
)
dev.off()

# Plot the prs_percentile  with colour by score
png("graphs/0001/0001_prevalence_percentile_colorbyscore.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
q <- ggplot(prs_percentile, aes(percentile, prevalence, color = (meanscore))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.001", 
       subtitle = "All samples, GCAT and Pandora \nColour by increase in score") +
  scale_fill_discrete(name="Mean Score") +
  theme_classic()
q + theme(
  plot.title = element_text(color = "dodgerblue4", size = 153 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic"),
  legend.title=element_blank())
dev.off()

# Plot the prs_percentile  with colour by score WITH LOG10
png("graphs/0001/0001_prevalence_percentile_colorbyscore_log10.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
a <- ggplot(prs_percentile, aes(percentile, prevalence, color = log10(meanscore))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.001", 
       subtitle = "All samples, GCAT and Pandora \nColour by increase in score") +
  scale_fill_discrete(name="Mean Score log10") +
  theme_classic()
a + theme(
  plot.title = element_text(color = "dodgerblue4", size = 153 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic"),
  legend.title=element_blank())
dev.off()

# Plot the prs_percentile  with colour by score WITH LOG10
png("graphs/0001/0001_prevalence_percentile_colorbythreshold.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
b <- ggplot(prs_percentile, aes(x = percentile, y = prevalence, color = as.factor(colorthreshold))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.001", 
       subtitle = "All samples, GCAT and Pandora \nColour by higher percentage of cases or control") +
  theme_classic()
b + scale_color_manual(values = c("lightskyblue", "dodgerblue4"),
                       name  ="",
                       breaks=c("0", "1"),
                       labels=c("High control \npercentage \n", "High case \npercentage")) +
  theme(
    plot.title = element_text(color = "dodgerblue4", size = 153 , face = "bold"),
    plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
    axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
    axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_hline(yintercept=50, linetype="dashed", color = "black") +
  geom_vline(xintercept=83.5, linetype="dashed", color = "black")
dev.off()



# Now we want to visualize GCAT and Pandora separately
head(general_score)
table(general_score$general)

general_score_gcat = general_score[which(general_score$general == 1),] 
general_score_pandora = general_score[which(general_score$general == 2),]


head(general_score_gcat)
dim(general_score_gcat)


# Create the prs_percentile variable
general_score_gcat$prs_percentile = 99
general_score_pandora$prs_percentile = 99

# Create the percentiles first number
my_seq_gcat = seq(1,nrow(general_score_gcat),round(nrow(general_score_gcat)*0.01))
my_seq_pandora = seq(1,nrow(general_score_pandora),round(nrow(general_score_pandora)*0.01))

# Fill the percentiles variable 
cont = 1
for(i in 1:(length(my_seq_gcat)-1)){
  
  general_score_gcat$prs_percentile[my_seq_gcat[i]:(my_seq_gcat[i+1]-1)] = cont
  
  cont = cont+1  
}
cont = 1
for(i in 1:(length(my_seq_pandora)-1)){
  
  general_score_pandora$prs_percentile[my_seq_pandora[i]:(my_seq_pandora[i+1]-1)] = cont
  
  cont = cont+1  
}

# Check the percentiles N
table(general_score_gcat$prs_percentile)
table(general_score_pandora$prs_percentile)
head(general_score_gcat)
head(general_score_pandora)

# Create the prs_percentiles dataframe 
prs_percentile_gcat = data.frame(percentile=seq(1,99,1))
prs_percentile_pandora = data.frame(percentile=seq(1,99,1))

# Create the prevalence variable
prs_percentile_gcat$meanscore = 0
prs_percentile_pandora$meanscore = 0

# Fill the prevalence GCAT and Pandora
for(i in 1:nrow(prs_percentile_gcat)){
  # Obtain the percentile information
  aux_G = general_score_gcat %>% filter(prs_percentile==i)
  aux_P = general_score_pandora %>% filter(prs_percentile==i)
  # Obtain the average score of the percentile
  prs_percentile_gcat$meanscore[i] = mean(aux_G$SCORE1_AVG)
  prs_percentile_pandora$meanscore[i] = mean(aux_P$SCORE1_AVG)
}

# Check 
prs_percentile_gcat
prs_percentile_pandora

# Change column names and merge the both datasets
colnames(prs_percentile_gcat)[2] = "meanscore_GCAT"
colnames(prs_percentile_pandora)[2] = "meanscore_Pandora"
prs_percentile_gcat_pandora = merge(prs_percentile_gcat, prs_percentile_pandora)
head(prs_percentile_gcat_pandora)


png("graphs/0001/0001_scoreincrease_bygroup.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
grob <- grobTree(textGrob("Score = 0.1 \npercentile = 14", x=0.007,  y=0.15, hjust=0,
                          gp=gpar(col="red", fontsize=7, fontface="italic")))
grob2 <- grobTree(textGrob("698 GCAT \n192 Pandora \n14%", x=0.01,  y=0.8, hjust=0,
                           gp=gpar(col="black", fontsize=8, fontface="italic")))
grob3 <- grobTree(textGrob("4290 GCAT \n1183 Pandora \n86%", x=0.2,  y=0.8, hjust=0,
                           gp=gpar(col="black", fontsize=8, fontface="italic")))
e <- ggplot(prs_percentile_gcat_pandora) +
  geom_line(aes(y= prs_percentile_gcat_pandora$meanscore_GCAT,
                x = prs_percentile_gcat_pandora$percentile),
            colour= "lightskyblue") +
  geom_line(aes(y= prs_percentile_gcat_pandora$meanscore_Pandora,
                x = prs_percentile_gcat_pandora$percentile),
            colour = "dodgerblue4") +
  labs(title = "Increase in Score by group PRS maf 0.001", 
       subtitle = "GCAT as ligthblue and Pandora as darkblue")+
  xlab("percentile")+
  ylab("PRS") +
  theme_classic()
e + theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_hline(yintercept=0.1, linetype="dashed", color = "red") +
  geom_vline(xintercept=14, linetype="dashed", color = "red") +
  annotation_custom(grob) + annotation_custom(grob2) + annotation_custom(grob3)
dev.off()



## Check the higher score GCAT participants conditions

# Read gcat
gcat_participants = readRDS("/imppc/labs/dnalab/share/PRS/GCAT_participants/gcat_participants_2019-03-06.rds")
head(gcat_participants)

# Read conditions
conditions = fread("/imppc/labs/dnalab/share/PRS/GCAT_participants/conditions.csv")
head(conditions)

# filter the antecedents 0, so the conditions will be of the patient not familiars
conditions = conditions %>% filter(antecedents==0)

# Merge the conditions and gcat_participants
conditions = left_join(conditions,gcat_participants %>% dplyr::select(entity_id,GENOTYPED_SAMPLE))

# Filter the GCAT participants
risk_score = general_score %>% filter(general == 1)

# Take the higher 50 score participants
risk_score = risk_score %>% arrange(desc(SCORE1_AVG)) %>% slice(1:50)

# Take the lower 50 score participants
#risk_score = risk_score %>% arrange(SCORE1_AVG) %>% slice(1:50)

# Filter the information from GCAT of this risk patients
gcat_risk = conditions[conditions$GENOTYPED_SAMPLE %in% risk_score$IID,]

head(gcat_risk)

table(gcat_risk$link)

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

