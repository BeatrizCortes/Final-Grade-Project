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
    plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
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
    plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
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
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
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
    plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
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

png("graphs/005/005_density.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
grob <- grobTree(textGrob("813 Pandora\n \n0 GCAT", x=0.7,  y=0.8, hjust=0,
                          gp=gpar(col="brown4", fontsize=7, fontface="italic")))
grob2 <- grobTree(textGrob("1058 \nPandora\n \n6 \nGCAT", x=0.55,  y=0.8, hjust=0,
                          gp=gpar(col="darkred", fontsize=7, fontface="italic")))
grob3 <- grobTree(textGrob("4960 \nPandora\n \n253 \nGCAT", x=0.4,  y=0.8, hjust=0,
                           gp=gpar(col="red2", fontsize=7, fontface="italic")))
grob4 <- grobTree(textGrob("169 \nPandora\n \n4896 \nGCAT", x=0.05,  y=0.8, hjust=0,
                           gp=gpar(col="brown1", fontsize=7, fontface="italic")))
e <- ggplot(general_score,aes(log10(SCORE1_AVG), fill = as.factor(general)))+
  geom_density(alpha = 0.2) +
  labs(title = "Increase in Score by group PRS maf 0.05")+
  xlab("Score log10")+
  ylab("") +
  theme_classic()
e + scale_fill_manual(values = c("cadetblue2", "blue4"),
                         name  ="",
                         breaks=c("1", "2"),
                         labels=c("GCAT \n", "Pandora")) +
  theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_vline(xintercept=general_score %>% 
               filter(general == 1) %>% 
               dplyr::select(SCORE1_AVG) %>% 
               max() %>% log10(), 
             linetype="dashed", color = "brown4") +
  geom_vline(xintercept= 0.5, linetype="dashed", color = "darkred") +
  geom_vline(xintercept= 0, linetype="dashed", color = "red2") +
  geom_vline(xintercept= -0.35, linetype="dashed", color = "brown1") +
  annotation_custom(grob) + annotation_custom(grob2) + 
  annotation_custom(grob3) + annotation_custom(grob4) 
dev.off()

# The maximun GCAT log10 score is 0.9404333
general_score %>% filter(general == 1) %>% dplyr::select(SCORE1_AVG) %>% max() %>% log10()

# Over log10(score) == 0.945 (max GCAT) only have Pandora
table(general_score[which(log10(general_score$SCORE1_AVG) > 0.945),]$general)

# Over 0.5 log10 score
table(general_score[which(log10(general_score$SCORE1_AVG) > 0.5),]$general)

# Lower than 0 log10 score
table(general_score[which(log10(general_score$SCORE1_AVG) < 0),]$general)

# Lower than 0 log10 score
table(general_score[which(log10(general_score$SCORE1_AVG) < -0.35),]$general)

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
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
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
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
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
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
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
    plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
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


## Histogram plot

png("graphs/001/001_density.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
grob <- grobTree(textGrob("883 Pandora\n \n0 GCAT", x=0.7,  y=0.8, hjust=0,
                          gp=gpar(col="brown4", fontsize=7, fontface="italic")))
grob2 <- grobTree(textGrob("1141 \nPandora\n \n11 \nGCAT", x=0.5,  y=0.8, hjust=0,
                           gp=gpar(col="darkred", fontsize=7, fontface="italic")))
grob3 <- grobTree(textGrob("128 \nPandora\n \n4904 \nGCAT", x=0.03,  y=0.8, hjust=0,
                           gp=gpar(col="red2", fontsize=7, fontface="italic")))
grob4 <- grobTree(textGrob("128 \nPandora\n \n4904 \nGCAT", x=0.03,  y=0.8, hjust=0,
                           gp=gpar(col="red2", fontsize=7, fontface="italic")))
e <- ggplot(general_score,aes(log10(SCORE1_AVG), fill = as.factor(general)))+
  geom_density(alpha = 0.2) +
  labs(title = "Increase in Score by group PRS maf 0.01")+
  xlab("Score log10")+
  ylab("") +
  theme_classic()
e + scale_fill_manual(values = c("cadetblue2", "blue4"),
                      name  ="",
                      breaks=c("1", "2"),
                      labels=c("GCAT \n", "Pandora")) +
  theme(
    plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
    plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
    axis.title.x = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_vline(xintercept= 0.54, 
             linetype="dashed", color = "brown4") +
  geom_vline(xintercept= 0, linetype="dashed", color = "darkred") +
  geom_vline(xintercept= -0.65, linetype="dashed", color = "red2") +
  annotation_custom(grob) + annotation_custom(grob2) + 
  annotation_custom(grob3) + annotation_custom(grob4)
dev.off()

# Over log10(score) == 0.54 (max GCAT) only have Pandora
table(general_score[which(log10(general_score$SCORE1_AVG) > 0.54),]$general)

# Lower than 0 log10 score
table(general_score[which(log10(general_score$SCORE1_AVG) > 0),]$general)

# Lower than -0.65 log10 score
table(general_score[which(log10(general_score$SCORE1_AVG) < -0.65),]$general)

table(general_score[which(log10(general_score$SCORE1_AVG) < 0 & log10(general_score$SCORE1_AVG) > -0.65),]$general)


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
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
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
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
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
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
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
    plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
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


## Histogram plot

png("graphs/0001/0001_density.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
grob <- grobTree(textGrob("1136 Pandora\n \n7 GCAT", x=0.5,  y=0.8, hjust=0,
                          gp=gpar(col="brown4", fontsize=7, fontface="italic")))
grob2 <- grobTree(textGrob("103 \nPandora\n \n39 \nGCAT", x=0.27,  y=0.8, hjust=0,
                           gp=gpar(col="black", fontsize=7, fontface="italic")))
grob4 <- grobTree(textGrob("136 \nPandora\n \n4942 \nGCAT", x=0.01,  y=0.8, hjust=0,
                           gp=gpar(col="brown1", fontsize=7, fontface="italic")))
e <- ggplot(general_score,aes(log10(SCORE1_AVG), fill = as.factor(general)))+
  geom_density(alpha = 0.2) +
  labs(title = "Increase in Score by group PRS maf 0.001")+
  xlab("Score log10")+
  ylab("") +
  theme_classic()
e + scale_fill_manual(values = c("cadetblue2", "blue4"),
                      name  ="",
                      breaks=c("1", "2"),
                      labels=c("GCAT \n", "Pandora")) +
  theme(
    plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
    plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
    axis.title.x = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_vline(xintercept=-0.7, 
             linetype="dashed", color = "brown4") +
  geom_vline(xintercept= -1.2, linetype="dashed", color = "brown1") +
  annotation_custom(grob) + annotation_custom(grob2) + annotation_custom(grob4)  
dev.off()

# Over -0.7 log10 score
table(general_score[which(log10(general_score$SCORE1_AVG) > -0.7),]$general)

# Lower than -1.2 log10 score
table(general_score[which(log10(general_score$SCORE1_AVG) < -1.2),]$general)

# Lower than 0 log10 score
table(general_score[which(log10(general_score$SCORE1_AVG) < -0.7 & log10(general_score$SCORE1_AVG) > -1.2),]$general)

ggplot(general_score,aes(log10(SCORE1_AVG), colour = as.factor(general)))+
  geom_freqpoly(binwidth = 0.05)


#############################################################################################
# BEST PRS (0.05) MODEL PREDICTION ON GROUPS                                                #
#############################################################################################
# p_value_threshold = 0.05, Pruning YES, Classification_threshold = 0.05

# Read the information of the different groups
psam = fread("../05_GWAS/outputs_plink/first_final_subset.psam")
# delete the columns we are not interested in
psam = psam %>% dplyr::select(IID,ovarian,breast,colon,pheno_30)

# Check
head(psam)
table(psam$ovarian) # 1 is control and 2 is case

# We want to obtain only the ones we have present in the 30%
psam = psam %>% filter(!is.na(pheno_30))
psam$pheno_30 = NULL

# Now we check: 1473 control in all groups
table(psam$ovarian) #135 cases
table(psam$breast) #253 cases
table(psam$colon) #114 cases

# Now we can calculate the score

# Create the variable to save the results
groups = c("ovarian", "breast", "colon")
prs_bygroup = data.frame(groups)
prs_bygroup

#Read 30 scores
general_score_30 = fread(paste0("outputs/prs_0.05_1e-06_YES_0.07_30.sscore"))
colnames(general_score_30)[3] = "general"
general_score_30 = general_score_30 %>% filter(!is.na(general))

# Merge the datasets 
score30_bygroup = merge(general_score_30, psam)
head(score30_bygroup)
head(psam)

#############################################################################################
# BEST PRS (0.05) ON OVARIAN                                                                #
#############################################################################################

# Select only the data for ovarian: control and cases 
score30_ovarian = score30_bygroup %>% filter(!is.na(ovarian)) %>% dplyr::select(IID, ovarian, NMISS_ALLELE_CT, SCORE1_AVG)
head(score30_ovarian)
dim(score30_ovarian)

# Now order by score
score30_ovarian = score30_ovarian %>% arrange(SCORE1_AVG)

# Create the prs_percentile variable
score30_ovarian$prs_percentile = 99

# Create the percentiles first number
my_seq = seq(1,nrow(score30_ovarian),round(nrow(score30_ovarian)*0.01))

# Fill the percentiles variable 
cont = 1
for(i in 1:(length(my_seq)-1)){
  
  score30_ovarian$prs_percentile[my_seq[i]:(my_seq[i+1]-1)] = cont
  
  cont = cont+1  
}

# Check the percentiles N
table(score30_ovarian$prs_percentile)

# Create the prs_percentiles dataframe 
prs_percentile = data.frame(percentile=seq(1,99,1))

# Create the prevalence variable
prs_percentile$prevalence = 0
prs_percentile$meanscore = 0

# Fill the prevalence
for(i in 1:nrow(prs_percentile)){
  # Obtain the percentile information
  aux = score30_ovarian %>% filter(prs_percentile==i)
  # Case and control information
  aux2 = prop.table(table(aux$ovarian))
  # Obtain the average score of the percentile
  prs_percentile$meanscore[i] = mean(aux$SCORE1_AVG)
  # Fill the prevalence
  if (dim(aux2) == 1){
    if (aux$ovarian[1]== 1){
      prs_percentile$prevalence[i] = 0
    }
    if (aux$ovarian[1]== 2){
      prs_percentile$prevalence[i] = 100
    }
  }
  if (dim(aux2) == 2){
    prs_percentile$prevalence[i] = aux2[2]*100 
  }
}

# Add information for colour
prs_percentile$colorthreshold = 0
prs_percentile[92:99,]$colorthreshold = 1


# Check 
prs_percentile

# Plot the prs_percentile 
png("graphs/bygroup/005_prevalence_percentile_ovarian.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
p <-  ggplot(prs_percentile, aes(percentile, prevalence)) +
  geom_point(color = "dodgerblue3") +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.05 \nOvarian", 
       subtitle = "All samples, GCAT and Pandora")+
  theme_classic() 
p + theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")
)
dev.off()

# Plot the prs_percentile  with colour by score
png("graphs/bygroup/005_prevalence_percentile_colorbyscore_ovarian.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
q <- ggplot(prs_percentile, aes(percentile, prevalence, color = (meanscore))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.05 \nOvarian", 
       subtitle = "All samples, GCAT and Pandora \nColour by increase in score") +
  scale_fill_discrete(name="Mean Score") +
  theme_classic()
q + theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic"),
  legend.title=element_blank())
dev.off()

# Plot the prs_percentile  with colour by score WITH LOG10
png("graphs/bygroup/005_prevalence_percentile_colorbyscore_log10_ovarian.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
a <- ggplot(prs_percentile, aes(percentile, prevalence, color = log10(meanscore))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.05 \nOvarian", 
       subtitle = "All samples, GCAT and Pandora \nColour by increase in score") +
  scale_fill_discrete(name="Mean Score log10") +
  theme_classic()
a + theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic"),
  legend.title=element_blank())
dev.off()

# Plot the prs_percentile  with colour by score WITH LOG10
png("graphs/bygroup/005_prevalence_percentile_colorbythreshold_ovarian.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
b <- ggplot(prs_percentile, aes(x = percentile, y = prevalence, color = as.factor(colorthreshold))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.05 \nOvarian", 
       subtitle = "All samples, GCAT and Pandora \nColour by higher percentage of cases or control") +
  theme_classic()
b + scale_color_manual(values = c("lightskyblue", "dodgerblue4"),
                       name  ="",
                       breaks=c("0", "1"),
                       labels=c("High control \npercentage \n", "High case \npercentage")) +
  theme(
    plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
    plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
    axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
    axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_hline(yintercept=18, linetype="dashed", color = "black") +
  geom_vline(xintercept=91.5, linetype="dashed", color = "black")
dev.off()

# Now we want to visualize GCAT and Pandora separately
head(score30_ovarian)
table(score30_ovarian$ovarian)

score30_ovarian_gcat = score30_ovarian[which(score30_ovarian$ovarian == 1),] 
score30_ovarian_pandora = score30_ovarian[which(score30_ovarian$ovarian == 2),]


head(score30_ovarian_gcat)
dim(score30_ovarian_gcat)


# Create the prs_percentile variable
score30_ovarian_gcat$prs_percentile = 99
score30_ovarian_pandora$prs_percentile = 99

# Create the percentiles first number
my_seq_gcat = seq(1,nrow(score30_ovarian_gcat),round(nrow(score30_ovarian_gcat)*0.01))
my_seq_pandora = seq(1,nrow(score30_ovarian_pandora),round(nrow(score30_ovarian_pandora)*0.01))

# Fill the percentiles variable 
cont = 1
for(i in 1:(length(my_seq_gcat)-1)){
  
  score30_ovarian_gcat$prs_percentile[my_seq_gcat[i]:(my_seq_gcat[i+1]-1)] = cont
  
  cont = cont+1  
}
cont = 1
for(i in 1:(length(my_seq_pandora)-1)){
  
  score30_ovarian_pandora$prs_percentile[my_seq_pandora[i]:(my_seq_pandora[i+1]-1)] = cont
  
  cont = cont+1  
}

# Check the percentiles N
table(score30_ovarian_gcat$prs_percentile)
table(score30_ovarian_pandora$prs_percentile)
head(score30_ovarian_gcat)
head(score30_ovarian_pandora)

# Create the prs_percentiles dataframe 
prs_percentile_gcat = data.frame(percentile=seq(1,99,1))
prs_percentile_pandora = data.frame(percentile=seq(1,99,1))

# Create the prevalence variable
prs_percentile_gcat$meanscore = 0
prs_percentile_pandora$meanscore = 0

# Fill the prevalence GCAT and Pandora
for(i in 1:nrow(prs_percentile_gcat)){
  # Obtain the percentile information
  aux_G = score30_ovarian_gcat %>% filter(prs_percentile==i)
  aux_P = score30_ovarian_pandora %>% filter(prs_percentile==i)
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


png("graphs/bygroup/005_scoreincrease_bygroup_ovarian.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
grob <- grobTree(textGrob("Score = 5 \npercentile = 33", x=0.007,  y=0.2, hjust=0,
                          gp=gpar(col="red", fontsize=7, fontface="italic")))
grob2 <- grobTree(textGrob("486 GCAT \n45 Pandora \n33%", x=0.15,  y=0.8, hjust=0,
                           gp=gpar(col="black", fontsize=8, fontface="italic")))
grob3 <- grobTree(textGrob("987 GCAT \n90 Pandora \n67%", x=0.4,  y=0.8, hjust=0,
                           gp=gpar(col="black", fontsize=8, fontface="italic")))
e <- ggplot(prs_percentile_gcat_pandora) +
  geom_line(aes(y= prs_percentile_gcat_pandora$meanscore_GCAT,
                x = prs_percentile_gcat_pandora$percentile),
            colour= "lightskyblue") +
  geom_line(aes(y= prs_percentile_gcat_pandora$meanscore_Pandora,
                x = prs_percentile_gcat_pandora$percentile),
            colour = "dodgerblue4") +
  labs(title = "Increase in Score by group PRS maf 0.05 \nOvarian", 
       subtitle = "GCAT as ligthblue and Pandora as darkblue")+
  xlab("percentile")+
  ylab("PRS") +
  theme_classic()
e + theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_hline(yintercept=5.1, linetype="dashed", color = "red") +
  geom_vline(xintercept=33, linetype="dashed", color = "red") +
  annotation_custom(grob) + annotation_custom(grob2) + annotation_custom(grob3)
dev.off()

head(max(prs_percentile_gcat_pandora$meanscore_GCAT))


## Histogram plot

png("graphs/bygroup/005_density_ovarian.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
grob <- grobTree(textGrob("109 Pandora (80%)\n \n2 GCAT", x=0.55,  y=0.8, hjust=0,
                          gp=gpar(col="brown4", fontsize=7, fontface="italic")))
grob2 <- grobTree(textGrob("11 Pandora\n22 GCAT", x=0.55,  y=0.40, hjust=0,
                           gp=gpar(col="black", fontsize=7, fontface="italic")))
grob3 <- grobTree(textGrob("116 \nPandora\n \n9 \nGCAT \n(86%)", x=0.358,  y=0.8, hjust=0,
                           gp=gpar(col="brown3", fontsize=7, fontface="italic")))
grob4 <- grobTree(textGrob("15 \nPandora\n \n1449 \nGCAT \n(98%)", x=0.02,  y=0.8, hjust=0,
                           gp=gpar(col="brown1", fontsize=7, fontface="italic")))
e <- ggplot(score30_ovarian,aes(log10(SCORE1_AVG), fill = as.factor(ovarian)))+
  geom_density(alpha = 0.2) +
  labs(title = "Increase in Score by group PRS maf 0.05 \nOvarian")+
  xlab("Score log10")+
  ylab("") +
  theme_classic()
e + scale_fill_manual(values = c("cadetblue2", "blue4"),
                      name  ="",
                      breaks=c("1", "2"),
                      labels=c("GCAT \n", "Pandora")) +
  theme(
    plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
    plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
    axis.title.x = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_vline(xintercept=-0.3, linetype="dashed", color = "brown1") +
  geom_vline(xintercept= -0.1, linetype="dashed", color = "brown3") +
  geom_vline(xintercept= 0.4, linetype="dashed", color = "brown4") +
  annotate("segment", x = -0.2, xend = 0.5, y = 1.55, yend = 1.55, colour = "black", 
           size=1, alpha=1, arrow=arrow())+
  annotation_custom(grob) + annotation_custom(grob2) + 
  annotation_custom(grob3) + annotation_custom(grob4) 
dev.off()

table(score30_ovarian[which(log10(score30_ovarian$SCORE1_AVG) > 0.4),]$ovarian) # 2 GCAT 109 Pandora
table(score30_ovarian[which(log10(score30_ovarian$SCORE1_AVG) < -0.3),]$ovarian) # 1449 GCAT 15 Pandora
table(score30_ovarian[which(log10(score30_ovarian$SCORE1_AVG) < -0.1 & log10(score30_ovarian$SCORE1_AVG) > -0.3),]$ovarian)#22 GCAT 11 Pandora
table(score30_ovarian[which(log10(score30_ovarian$SCORE1_AVG) > -0.1),]$ovarian) # 9 GCAT 116 Pandora


#############################################################################################
# BEST PRS (0.05) ON BREAST                                                                 #
#############################################################################################

# Select only the data for breast: control and cases 
score30_breast = score30_bygroup %>% filter(!is.na(breast)) %>% dplyr::select(IID, breast, NMISS_ALLELE_CT, SCORE1_AVG)
head(score30_breast)
dim(score30_breast)

# Now order by score
score30_breast = score30_breast %>% arrange(SCORE1_AVG)

# Create the prs_percentile variable
score30_breast$prs_percentile = 99

# Create the percentiles first number
my_seq = seq(1,nrow(score30_breast),round(nrow(score30_breast)*0.01))

# Fill the percentiles variable 
cont = 1
for(i in 1:(length(my_seq)-1)){
  
  score30_breast$prs_percentile[my_seq[i]:(my_seq[i+1]-1)] = cont
  
  cont = cont+1  
}

# Check the percentiles N
table(score30_breast$prs_percentile)

# Create the prs_percentiles dataframe 
prs_percentile = data.frame(percentile=seq(1,99,1))

# Create the prevalence variable
prs_percentile$prevalence = 0
prs_percentile$meanscore = 0

# Fill the prevalence
for(i in 1:nrow(prs_percentile)){
  # Obtain the percentile information
  aux = score30_breast %>% filter(prs_percentile==i)
  # Case and control information
  aux2 = prop.table(table(aux$breast))
  # Obtain the average score of the percentile
  prs_percentile$meanscore[i] = mean(aux$SCORE1_AVG)
  # Fill the prevalence
  if (dim(aux2) == 1){
    if (aux$breast[1]== 1){
      prs_percentile$prevalence[i] = 0
    }
    if (aux$breast[1]== 2){
      prs_percentile$prevalence[i] = 100
    }
  }
  if (dim(aux2) == 2){
    prs_percentile$prevalence[i] = aux2[2]*100 
  }
}

# Add information for colour
prs_percentile$colorthreshold = 0
prs_percentile[87:99,]$colorthreshold = 1


# Check 
prs_percentile

# Plot the prs_percentile 
png("graphs/bygroup/005_prevalence_percentile_breast.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
p <-  ggplot(prs_percentile, aes(percentile, prevalence)) +
  geom_point(color = "dodgerblue3") +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.05 \nBreast", 
       subtitle = "All samples, GCAT and Pandora")+
  theme_classic() 
p + theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic"))
dev.off()

# Plot the prs_percentile  with colour by score
png("graphs/bygroup/005_prevalence_percentile_colorbyscore_breast.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
q <- ggplot(prs_percentile, aes(percentile, prevalence, color = (meanscore))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.05 \nBreast", 
       subtitle = "All samples, GCAT and Pandora \nColour by increase in score") +
  scale_fill_discrete(name="Mean Score") +
  theme_classic()
q + theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic"),
  legend.title=element_blank())
dev.off()

# Plot the prs_percentile  with colour by score WITH LOG10
png("graphs/bygroup/005_prevalence_percentile_colorbyscore_log10_breast.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
a <- ggplot(prs_percentile, aes(percentile, prevalence, color = log10(meanscore))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.05 \nBreast", 
       subtitle = "All samples, GCAT and Pandora \nColour by increase in score") +
  scale_fill_discrete(name="Mean Score log10") +
  theme_classic()
a + theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic"),
  legend.title=element_blank())
dev.off()

# Plot the prs_percentile  with colour by score WITH LOG10
png("graphs/bygroup/005_prevalence_percentile_colorbythreshold_breast.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
b <- ggplot(prs_percentile, aes(x = percentile, y = prevalence, color = as.factor(colorthreshold))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.05 \nBreast", 
       subtitle = "All samples, GCAT and Pandora \nColour by higher percentage of cases or control") +
  theme_classic()
b + scale_color_manual(values = c("lightskyblue", "dodgerblue4"),
                       name  ="",
                       breaks=c("0", "1"),
                       labels=c("High control \npercentage \n", "High case \npercentage")) +
  theme(
    plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
    plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
    axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
    axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_hline(yintercept=25, linetype="dashed", color = "black") +
  geom_vline(xintercept=86.5, linetype="dashed", color = "black")
dev.off()

# Now we want to visualize GCAT and Pandora separately
head(score30_breast)
table(score30_breast$breast)

score30_breast_gcat = score30_breast[which(score30_breast$breast == 1),] 
score30_breast_pandora = score30_breast[which(score30_breast$breast == 2),]


head(score30_breast_gcat)
dim(score30_breast_gcat)


# Create the prs_percentile variable
score30_breast_gcat$prs_percentile = 99
score30_breast_pandora$prs_percentile = 99

# Create the percentiles first number
my_seq_gcat = seq(1,nrow(score30_breast_gcat),round(nrow(score30_breast_gcat)*0.01))
my_seq_pandora = seq(1,nrow(score30_breast_pandora),round(nrow(score30_breast_pandora)*0.01))

# Fill the percentiles variable 
cont = 1
for(i in 1:(length(my_seq_gcat)-1)){
  
  score30_breast_gcat$prs_percentile[my_seq_gcat[i]:(my_seq_gcat[i+1]-1)] = cont
  
  cont = cont+1  
}
cont = 1
for(i in 1:(length(my_seq_pandora)-1)){
  
  score30_breast_pandora$prs_percentile[my_seq_pandora[i]:(my_seq_pandora[i+1]-1)] = cont
  
  cont = cont+1  
}

# Check the percentiles N
table(score30_breast_gcat$prs_percentile)
table(score30_breast_pandora$prs_percentile)
head(score30_breast_gcat)
head(score30_breast_pandora)

# Create the prs_percentiles dataframe 
prs_percentile_gcat = data.frame(percentile=seq(1,99,1))
prs_percentile_pandora = data.frame(percentile=seq(1,99,1))

# Create the prevalence variable
prs_percentile_gcat$meanscore = 0
prs_percentile_pandora$meanscore = 0

# Fill the prevalence GCAT and Pandora
for(i in 1:nrow(prs_percentile_gcat)){
  # Obtain the percentile information
  aux_G = score30_breast_gcat %>% filter(prs_percentile==i)
  aux_P = score30_breast_pandora %>% filter(prs_percentile==i)
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


png("graphs/bygroup/005_scoreincrease_bygroup_breast.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
grob <- grobTree(textGrob("Score = 5 \npercentile = 21", x=0.007,  y=0.2, hjust=0,
                          gp=gpar(col="red", fontsize=7, fontface="italic")))
grob2 <- grobTree(textGrob("1164 GCAT \n53 Pandora \n21%", x=0.05,  y=0.8, hjust=0,
                           gp=gpar(col="black", fontsize=8, fontface="italic")))
grob3 <- grobTree(textGrob("309 GCAT \n200 Pandora \n79%", x=0.28,  y=0.8, hjust=0,
                           gp=gpar(col="black", fontsize=8, fontface="italic")))
e <- ggplot(prs_percentile_gcat_pandora) +
  geom_line(aes(y= prs_percentile_gcat_pandora$meanscore_GCAT,
                x = prs_percentile_gcat_pandora$percentile),
            colour= "lightskyblue") +
  geom_line(aes(y= prs_percentile_gcat_pandora$meanscore_Pandora,
                x = prs_percentile_gcat_pandora$percentile),
            colour = "dodgerblue4") +
  labs(title = "Increase in Score by group PRS maf 0.05 \nBreast", 
       subtitle = "GCAT as ligthblue and Pandora as darkblue")+
  xlab("percentile")+
  ylab("PRS") +
  theme_classic()
e + theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_hline(yintercept=5.1, linetype="dashed", color = "red") +
  geom_vline(xintercept=21, linetype="dashed", color = "red") +
  annotation_custom(grob) + annotation_custom(grob2) + annotation_custom(grob3)
dev.off()

head(max(prs_percentile_gcat_pandora$meanscore_GCAT))


## Histogram plot

png("graphs/bygroup/005_density_breast.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
grob <- grobTree(textGrob("200 Pandora (79%)\n \n2 GCAT", x=0.55,  y=0.8, hjust=0,
                          gp=gpar(col="brown4", fontsize=7, fontface="italic")))
grob2 <- grobTree(textGrob("6 Pandora\n15 GCAT", x=0.55,  y=0.40, hjust=0,
                           gp=gpar(col="black", fontsize=7, fontface="italic")))
grob3 <- grobTree(textGrob("216 \nPandora\n \n9 \nGCAT \n(85%)", x=0.358,  y=0.8, hjust=0,
                           gp=gpar(col="brown3", fontsize=7, fontface="italic")))
grob4 <- grobTree(textGrob("31 \nPandora\n \n1449 \nGCAT \n(98%)", x=0.02,  y=0.8, hjust=0,
                           gp=gpar(col="brown1", fontsize=7, fontface="italic")))

e <- ggplot(score30_breast,aes(log10(SCORE1_AVG), fill = as.factor(breast)))+
  geom_density(alpha = 0.2) +
  labs(title = "Increase in Score by group PRS maf 0.05 \nBreast")+
  xlab("Score log10")+
  ylab("") +
  theme_classic()
e + scale_fill_manual(values = c("cadetblue2", "blue4"),
                      name  ="",
                      breaks=c("1", "2"),
                      labels=c("GCAT \n", "Pandora")) +
  theme(
    plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
    plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
    axis.title.x = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_vline(xintercept=-0.3, linetype="dashed", color = "brown1") +
  geom_vline(xintercept= -0.1, linetype="dashed", color = "brown3") +
  geom_vline(xintercept= 0.4, linetype="dashed", color = "brown4") +
  annotate("segment", x = -0.2, xend = 0.5, y = 1.55, yend = 1.55, colour = "black", 
           size=1, alpha=1, arrow=arrow())+
  annotation_custom(grob) + annotation_custom(grob2) + 
  annotation_custom(grob3) + annotation_custom(grob4) 
dev.off()

table(score30_breast[which(log10(score30_breast$SCORE1_AVG) > 0.4),]$breast) # 2 GCAT 200 Pandora
table(score30_breast[which(log10(score30_breast$SCORE1_AVG) < -0.3),]$breast) # 1449 GCAT 31 Pandora
table(score30_breast[which(log10(score30_breast$SCORE1_AVG) < -0.1 & log10(score30_breast$SCORE1_AVG) > -0.3),]$breast)#15 GCAT 6 Pandora
table(score30_breast[which(log10(score30_breast$SCORE1_AVG) > -0.1),]$breast) # 9 GCAT 216 Pandora


#############################################################################################
# BEST PRS (0.05) ON COLON                                                                 #
#############################################################################################

# Select only the data for colon: control and cases 
score30_colon = score30_bygroup %>% filter(!is.na(colon)) %>% dplyr::select(IID, colon, NMISS_ALLELE_CT, SCORE1_AVG)
head(score30_colon)
dim(score30_colon)

# Now order by score
score30_colon = score30_colon %>% arrange(SCORE1_AVG)

# Create the prs_percentile variable
score30_colon$prs_percentile = 99

# Create the percentiles first number
my_seq = seq(1,nrow(score30_colon),round(nrow(score30_colon)*0.01))

# Fill the percentiles variable 
cont = 1
for(i in 1:(length(my_seq)-1)){
  
  score30_colon$prs_percentile[my_seq[i]:(my_seq[i+1]-1)] = cont
  
  cont = cont+1  
}

# Check the percentiles N
table(score30_colon$prs_percentile)

# Create the prs_percentiles dataframe 
prs_percentile = data.frame(percentile=seq(1,99,1))

# Create the prevalence variable
prs_percentile$prevalence = 0
prs_percentile$meanscore = 0

# Fill the prevalence
for(i in 1:nrow(prs_percentile)){
  # Obtain the percentile information
  aux = score30_colon %>% filter(prs_percentile==i)
  # Case and control information
  aux2 = prop.table(table(aux$colon))
  # Obtain the average score of the percentile
  prs_percentile$meanscore[i] = mean(aux$SCORE1_AVG)
  # Fill the prevalence
  if (dim(aux2) == 1){
    if (aux$colon[1]== 1){
      prs_percentile$prevalence[i] = 0
    }
    if (aux$colon[1]== 2){
      prs_percentile$prevalence[i] = 100
    }
  }
  if (dim(aux2) == 2){
    prs_percentile$prevalence[i] = aux2[2]*100 
  }
}

# Add information for colour
prs_percentile$colorthreshold = 0
prs_percentile[93:99,]$colorthreshold = 1


# Check 
prs_percentile

# Plot the prs_percentile 
png("graphs/bygroup/005_prevalence_percentile_colon.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
p <-  ggplot(prs_percentile, aes(percentile, prevalence)) +
  geom_point(color = "dodgerblue3") +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.05 \nColon", 
       subtitle = "All samples, GCAT and Pandora")+
  theme_classic() 
p + theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic"))
dev.off()

# Plot the prs_percentile  with colour by score
png("graphs/bygroup/005_prevalence_percentile_colorbyscore_colon.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
q <- ggplot(prs_percentile, aes(percentile, prevalence, color = (meanscore))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.05 \nColon", 
       subtitle = "All samples, GCAT and Pandora \nColour by increase in score") +
  scale_fill_discrete(name="Mean Score") +
  theme_classic()
q + theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic"),
  legend.title=element_blank())
dev.off()

# Plot the prs_percentile  with colour by score WITH LOG10
png("graphs/bygroup/005_prevalence_percentile_colorbyscore_log10_colon.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
a <- ggplot(prs_percentile, aes(percentile, prevalence, color = log10(meanscore))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.05 \nColon", 
       subtitle = "All samples, GCAT and Pandora \nColour by increase in score") +
  scale_fill_discrete(name="Mean Score log10") +
  theme_classic()
a + theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic"),
  legend.title=element_blank())
dev.off()

# Plot the prs_percentile  with colour by score WITH LOG10
png("graphs/bygroup/005_prevalence_percentile_colorbythreshold_colon.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
b <- ggplot(prs_percentile, aes(x = percentile, y = prevalence, color = as.factor(colorthreshold))) +
  geom_point() +
  labs(title = "Prevalence of cancer by percentile PRS maf 0.05 \nColon", 
       subtitle = "All samples, GCAT and Pandora \nColour by higher percentage of cases or control") +
  theme_classic()
b + scale_color_manual(values = c("lightskyblue", "dodgerblue4"),
                       name  ="",
                       breaks=c("0", "1"),
                       labels=c("High control \npercentage \n", "High case \npercentage")) +
  theme(
    plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
    plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
    axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
    axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_hline(yintercept=50, linetype="dashed", color = "black") +
  geom_vline(xintercept=92.5, linetype="dashed", color = "black")
dev.off()

# Now we want to visualize GCAT and Pandora separately
head(score30_colon)
table(score30_colon$colon)

score30_colon_gcat = score30_colon[which(score30_colon$colon == 1),] 
score30_colon_pandora = score30_colon[which(score30_colon$colon == 2),]


head(score30_colon_gcat)
dim(score30_colon_gcat)


# Create the prs_percentile variable
score30_colon_gcat$prs_percentile = 99
score30_colon_pandora$prs_percentile = 99

# Create the percentiles first number
my_seq_gcat = seq(1,nrow(score30_colon_gcat),round(nrow(score30_colon_gcat)*0.01))
my_seq_pandora = seq(1,nrow(score30_colon_pandora),round(nrow(score30_colon_pandora)*0.01))

# Fill the percentiles variable 
cont = 1
for(i in 1:(length(my_seq_gcat)-1)){
  
  score30_colon_gcat$prs_percentile[my_seq_gcat[i]:(my_seq_gcat[i+1]-1)] = cont
  
  cont = cont+1  
}
cont = 1
for(i in 1:(length(my_seq_pandora)-1)){
  
  score30_colon_pandora$prs_percentile[my_seq_pandora[i]:(my_seq_pandora[i+1]-1)] = cont
  
  cont = cont+1  
}

# Check the percentiles N
table(score30_colon_gcat$prs_percentile)
table(score30_colon_pandora$prs_percentile)
head(score30_colon_gcat)
head(score30_colon_pandora)

# Create the prs_percentiles dataframe 
prs_percentile_gcat = data.frame(percentile=seq(1,99,1))
prs_percentile_pandora = data.frame(percentile=seq(1,99,1))

# Create the prevalence variable
prs_percentile_gcat$meanscore = 0
prs_percentile_pandora$meanscore = 0

# Fill the prevalence GCAT and Pandora
for(i in 1:nrow(prs_percentile_gcat)){
  # Obtain the percentile information
  aux_G = score30_colon_gcat %>% filter(prs_percentile==i)
  aux_P = score30_colon_pandora %>% filter(prs_percentile==i)
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


png("graphs/bygroup/005_scoreincrease_bygroup_colon.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
grob <- grobTree(textGrob("Score = 5 \npercentile = 21", x=0.007,  y=0.2, hjust=0,
                          gp=gpar(col="red", fontsize=7, fontface="italic")))
grob2 <- grobTree(textGrob("1164 GCAT \n24 Pandora \n21%", x=0.05,  y=0.8, hjust=0,
                           gp=gpar(col="black", fontsize=8, fontface="italic")))
grob3 <- grobTree(textGrob("309 GCAT \n90 Pandora \n79%", x=0.28,  y=0.8, hjust=0,
                           gp=gpar(col="black", fontsize=8, fontface="italic")))
e <- ggplot(prs_percentile_gcat_pandora) +
  geom_line(aes(y= prs_percentile_gcat_pandora$meanscore_GCAT,
                x = prs_percentile_gcat_pandora$percentile),
            colour= "lightskyblue") +
  geom_line(aes(y= prs_percentile_gcat_pandora$meanscore_Pandora,
                x = prs_percentile_gcat_pandora$percentile),
            colour = "dodgerblue4") +
  labs(title = "Increase in Score by group PRS maf 0.05 \nColon", 
       subtitle = "GCAT as ligthblue and Pandora as darkblue")+
  xlab("percentile")+
  ylab("PRS") +
  theme_classic()
e + theme(
  plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
  plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
  axis.title.x = element_text(color="dodgerblue4", size=10, face="italic"),
  axis.title.y = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_hline(yintercept=5.1, linetype="dashed", color = "red") +
  geom_vline(xintercept=21, linetype="dashed", color = "red") +
  annotation_custom(grob) + annotation_custom(grob2) + annotation_custom(grob3)
dev.off()

head(max(prs_percentile_gcat_pandora$meanscore_GCAT))


## Histogram plot

png("graphs/bygroup/005_density_colon.png", width = 1500, height = 1500, units = "px", res = 300)
grid.newpage()
grob <- grobTree(textGrob("85 Pandora (79%)\n \n2 GCAT", x=0.55,  y=0.8, hjust=0,
                          gp=gpar(col="brown4", fontsize=7, fontface="italic")))
grob2 <- grobTree(textGrob("6 Pandora\n15 GCAT", x=0.55,  y=0.40, hjust=0,
                           gp=gpar(col="black", fontsize=7, fontface="italic")))
grob3 <- grobTree(textGrob("96 \nPandora\n \n9 \nGCAT \n(85%)", x=0.358,  y=0.8, hjust=0,
                           gp=gpar(col="brown3", fontsize=7, fontface="italic")))
grob4 <- grobTree(textGrob("12 \nPandora\n \n1449 \nGCAT \n(98%)", x=0.02,  y=0.8, hjust=0,
                           gp=gpar(col="brown1", fontsize=7, fontface="italic")))

e <- ggplot(score30_colon,aes(log10(SCORE1_AVG), fill = as.factor(colon)))+
  geom_density(alpha = 0.2) +
  labs(title = "Increase in Score by group PRS maf 0.05 \nColon")+
  xlab("Score log10")+
  ylab("") +
  theme_classic()
e + scale_fill_manual(values = c("cadetblue2", "blue4"),
                      name  ="",
                      breaks=c("1", "2"),
                      labels=c("GCAT \n", "Pandora")) +
  theme(
    plot.title = element_text(color = "dodgerblue4", size = 13 , face = "bold"),
    plot.subtitle = element_text(color = "dodgerblue3", size = 12 , face = "italic"),
    axis.title.x = element_text(color="dodgerblue4", size=10, face="italic")) +
  geom_vline(xintercept=-0.3, linetype="dashed", color = "brown1") +
  geom_vline(xintercept= -0.1, linetype="dashed", color = "brown3") +
  geom_vline(xintercept= 0.4, linetype="dashed", color = "brown4") +
  annotate("segment", x = -0.2, xend = 0.5, y = 1.55, yend = 1.55, colour = "black", 
           size=1, alpha=1, arrow=arrow())+
  annotation_custom(grob) + annotation_custom(grob2) + 
  annotation_custom(grob3) + annotation_custom(grob4) 
dev.off()

table(score30_colon[which(log10(score30_colon$SCORE1_AVG) > 0.4),]$colon) # 2 GCAT 85 Pandora
table(score30_colon[which(log10(score30_colon$SCORE1_AVG) < -0.3),]$colon) # 1449 GCAT 12 Pandora
table(score30_colon[which(log10(score30_colon$SCORE1_AVG) < -0.1 & log10(score30_colon$SCORE1_AVG) > -0.3),]$colon)#15 GCAT 6 Pandora
table(score30_colon[which(log10(score30_colon$SCORE1_AVG) > -0.1),]$colon) # 9 GCAT 96 Pandora


