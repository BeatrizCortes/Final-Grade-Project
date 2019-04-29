######################################################################################
# 06 GWAS                                                                            #
######################################################################################

# Libraries
library(dplyr)
library(data.table)
#install.packages("qqman")
library(qqman)

# Set the working directory 
setwd("/imppc/labs/dnalab/share/PRS/06_GWAS/")

# Read the data 
psam = fread("../05_GWAS/outputs_plink/first_final_subset.psam")
pvar = fread("../05_GWAS/outputs_plink/first_final_subset.pvar")

# Check
dim(pvar) # 13023     5
dim(psam) # 6569      4

## ADD CASE AND CONTROL ##

# First we create the variable and put all as 1 (control)
psam$pheno = 1 

# Check
table(psam$pheno)

# But now we have to take into acount that inside pandora there are some GCAT controls 
# We need to change them again to 1 (case)

samples_information = fread("/imppc/labs/dnalab/share/PRS/02_Pandora_data_preparation/samples_information.txt")
samples_information_sel = fread("/imppc/labs/dnalab/share/PRS/02_Pandora_data_preparation/selected_merge_samples.txt")

dim(samples_information_sel) # 35893  1587 contain the 1581 selected samples plus 6 variables

head(samples_information) # We are interested into Type (Case/Control)

dim(samples_information) # 1946
(nrow(psam) - 4988) # 1581

# Now we take the ids of the selected samples
sel_samples = colnames(samples_information_sel[,7:ncol(samples_information_sel)])
head(sel_samples) # Check

# Now we use this to filter the samples_information to keep only the selected samples
samples_information = samples_information[samples_information$sampleId %in% sel_samples,]
dim(samples_information) # 1587    7
head(samples_information) # This samples have the same order as the genotypes file we use to create our data

# We have duplicated samplesId
length(unique(samples_information$sampleId))

# Check which:
samples_information[which(duplicated(samples_information$sampleId)),] # 5041

# Check the first sample ID 
samples_information[samples_information$sampleId == "5042"]

# We have to delete this duplicated 
samples_information = samples_information[!which(duplicated(samples_information$sampleId)),]

# Create the vector to sustitute 
table(samples_information$Type)
# Case GCAT 
#1375  206 
pheno_Pandora = samples_information$Type

# The Cases to 2 and Cases to 1 
pheno_Pandora = ifelse(pheno_Pandora == "Case", 2, 1)

# Check
table(pheno_Pandora)
#    1    2 
#  206 1375

# Now add it to psam
psam$pheno[4989:nrow(psam)] = pheno_Pandora

# Check the final result
table(psam$pheno)

## ADD CANCER TYPE ##
# We will add onw column for each cancer type and then indicate the value for each sample
# We have the following cancer types 
table(samples_information$cascades.polides)
# Control                -> Control 
# Breast                 -> Breast
# Polyposis              -> Colon
# BC/OC/HBOC + Other     -> Breast + Ovary 
# HNPCC                  -> Colon
# HBOC                   -> Breast + Ovary
# Ovary                  -> Ovary

# OVARIAN
# Create all as control
psam$ovarian = 1 

# Create the auxiliar variable
pheno_ovarian = samples_information$cascades.polides

# The Cases to 2 and Cases to 1 
pheno_ovarian = ifelse(pheno_ovarian == "BC/OC/HBOC + Other" | 
                         pheno_ovarian == "HBOC" |
                         pheno_ovarian == "Ovary", 2, NA)

table(pheno_ovarian) # 41+158+225 = 424 cases

# Now add it to psam
psam$ovarian[4989:nrow(psam)] = pheno_ovarian

# Check 
table(psam$ovarian)

# BREAST
# Create all as control
psam$breast = 1

# Create the auxiliar variable
pheno_breast = samples_information$cascades.polides

# The Cases to 2 and Cases to 1
pheno_breast = ifelse(pheno_breast == "BC/OC/HBOC + Other" | 
                         pheno_breast == "HBOC" |
                         pheno_breast == "Breast", 2, NA)

table(pheno_breast) # 41+158+632 = 831 cases

# Now add it to psam
psam$breast[4989:nrow(psam)] = pheno_breast

# Check 
table(psam$breast)

# COLON
# Create all as control
psam$colon = 1 

# Create the auxiliar variable
pheno_colon = samples_information$cascades.polides

# The Cases to 2 and Cases to 1 
pheno_colon = ifelse(pheno_colon == "HNPCC" | 
                        pheno_colon == "Polyposis", 2, NA)

table(pheno_colon) # 168+151 = 319 cases

# Now add it to psam
psam$colon[4989:nrow(psam)] = pheno_colon

# Check 
table(psam$colon)

# Write the psam 
write.table(psam,"../05_GWAS/outputs_plink/first_final_subset.psam",
            quote = F,row.names = F)

# GENERAL ANALYSIS # 
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --pheno-name pheno ",
              " --glm firth ", # help to order the betas to take into accound not balanced cases
              " --maf 0.01 ", 
              " --out /imppc/labs/dnalab/share/PRS/06_GWAS/outputs_plink/first_final_subset"))

gwas = fread("/imppc/labs/dnalab/share/PRS/06_GWAS/outputs_plink/first_final_subset.pheno.glm.firth")

head(gwas)

# Create the gwas_manhattan with the columns we need 
gwas_manhattan = gwas %>% dplyr::select(ID,`#CHROM`,POS,P)
# Change the columns name
colnames(gwas_manhattan) = c("SNP","CHR","BP","P")

# Make the general GWAS
# dir.create("graphs/General_GWAS")

png("graphs/General_GWAS/general.png")
  manhattan(gwas_manhattan)
  title("All samples")
dev.off()

# Loop to make the plots

# First we need to create the vector with the chromosomes names:
chromosomes = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                "22")

for (c in chromosomes){
  png(paste0("graphs/General_GWAS/Chr", c ,"_manhattan.png"))
  manhattan(gwas_manhattan %>% filter(CHR == as.integer(c)))
  dev.off()
}

# check the correlated
general_variants = gwas %>% filter(P<5e-8) #412

# OVARIAN ANALYSIS # 
# we want now analyze by OVARIAN
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --pheno-name ovarian ",
              " --glm firth ", # help to order the betas to take into accound not balanced cases
              " --maf 0.01 ", 
              " --out /imppc/labs/dnalab/share/PRS/06_GWAS/outputs_plink/first_final_subset"))

gwas_O = fread("/imppc/labs/dnalab/share/PRS/06_GWAS/outputs_plink/first_final_subset.ovarian.glm.firth")

head(gwas_O)

# Create the gwas_manhattan with the columns we need 
gwas_manhattan_O = gwas_O %>% dplyr::select(ID,`#CHROM`,POS,P)
# Change the columns name
colnames(gwas_manhattan_O) = c("SNP","CHR","BP","P")

# Make the general GWAS
# dir.create("graphs/Ovarian_GWAS")

png("graphs/Ovarian_GWAS/general_ovarian.png")
  manhattan(gwas_manhattan_O)
  title("All samples Ovarian")
dev.off()

# Loop to make the plots

for (c in chromosomes){
  png(paste0("graphs/Ovarian_GWAS/Chr", c ,"_manhattan_ovarian.png"))
  manhattan(gwas_manhattan_O %>% filter(CHR == as.integer(c)))
  dev.off()
}

# check the correlated
ovarian_variants = gwas_O %>% filter(P<5e-8) #188


# BREAST ANALYSIS # 
# we want now analyze by Breast
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --pheno-name breast ",
              " --glm firth ", # help to order the betas to take into accound not balanced cases
              " --maf 0.01 ", 
              " --out /imppc/labs/dnalab/share/PRS/06_GWAS/outputs_plink/first_final_subset"))

gwas_B = fread("/imppc/labs/dnalab/share/PRS/06_GWAS/outputs_plink/first_final_subset.breast.glm.firth")

head(gwas_B)

# Create the gwas_manhattan with the columns we need 
gwas_manhattan_C = gwas_B %>% dplyr::select(ID,`#CHROM`,POS,P)
# Change the columns name
colnames(gwas_manhattan_C) = c("SNP","CHR","BP","P")

# Make the general GWAS
# dir.create("graphs/Breast_GWAS")

png("graphs/Breast_GWAS/general_breast.png")
manhattan(gwas_manhattan_C)
title("All samples")
dev.off()

# Loop to make the plots

for (c in chromosomes){
  png(paste0("graphs/Breast_GWAS/Chr", c ,"_manhattan_breast.png"))
  manhattan(gwas_manhattan_C %>% filter(CHR == as.integer(c)))
  dev.off()
}

# check the correlated
breast_variants = gwas_B %>% filter(P<5e-8) #318


# COLON ANALYSIS # 
# we want now analyze by Colon
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --pheno-name colon ",
              " --glm firth ", # help to order the betas to take into accound not balanced cases
              " --maf 0.01 ", 
              " --out /imppc/labs/dnalab/share/PRS/06_GWAS/outputs_plink/first_final_subset"))

gwas_C = fread("/imppc/labs/dnalab/share/PRS/06_GWAS/outputs_plink/first_final_subset.colon.glm.firth")

head(gwas_C)

# Create the gwas_manhattan with the columns we need 
gwas_manhattan_C = gwas_C %>% dplyr::select(ID,`#CHROM`,POS,P)
# Change the columns name
colnames(gwas_manhattan_C) = c("SNP","CHR","BP","P")

# Make the general GWAS
# dir.create("graphs/Colon_GWAS")

png("graphs/Colon_GWAS/general_colon.png")
manhattan(gwas_manhattan_C)
title("All samples")
dev.off()

# Loop to make the plots

for (c in chromosomes){
  png(paste0("graphs/Colon_GWAS/Chr", c ,"_manhattan_colon.png"))
  manhattan(gwas_manhattan_C %>% filter(CHR == as.integer(c)))
  dev.off()
}

# check the correlated
colon_variants = gwas_C %>% filter(P<5e-8)#153


## Search the information of this selected SNPs ##

genes_information = fread("/imppc/labs/dnalab/share/PRS/02_Pandora_data_preparation/panels.txt")
genes_information = genes_information[,1:2]
head(genes_information)

# Create a common variable chr_pos
general_variants$chr_pos = paste0(general_variants$`#CHROM`, "_", general_variants$POS)
ovarian_variants$chr_pos = paste0(ovarian_variants$`#CHROM`, "_", ovarian_variants$POS)
breast_variants$chr_pos = paste0(breast_variants$`#CHROM`, "_", breast_variants$POS)
colon_variants$chr_pos = paste0(colon_variants$`#CHROM`, "_", colon_variants$POS)

# Check

length(general_variants$chr_pos) - length(which(general_variants$chr_pos %in% genes_information$chr_pos)) ### 2 IVAN
which(!general_variants$chr_pos %in% genes_information$chr_pos) # 117 156

length(ovarian_variants$chr_pos) - length(which(ovarian_variants$chr_pos %in% ovarian_variants$chr_pos))
which(!ovarian_variants$chr_pos %in% genes_information$chr_pos) # 76

length(breast_variants$chr_pos) - length(which(breast_variants$chr_pos %in% breast_variants$chr_pos))
which(!breast_variants$chr_pos %in% genes_information$chr_pos) # 94 123

length(colon_variants$chr_pos) - length(which(colon_variants$chr_pos %in% colon_variants$chr_pos))
which(!colon_variants$chr_pos %in% genes_information$chr_pos) # 49


# Create the gene variable
general_variants$Gene = NA
ovarian_variants$Gene = NA
breast_variants$Gene = NA
colon_variants$Gene = NA

## GENERAL GENE SEARCHING ##

sel = which(general_variants$chr_pos %in% genes_information$chr_pos)# Select the common
for (i in sel){
  print(i)
  chrpos = general_variants[i,]$chr_pos
  general_variants[i,]$Gene = unique(genes_information[which(genes_information$chr_pos == chrpos),]$Gene)
}

# Check and manual fill GENERAL VARIANTS
general_variants[which(is.na(general_variants$Gene)),]$chr_pos # "5_60240986"  "7_142457132"

# Search gene
gene = unique(genes_information[which(genes_information$chr_pos > "5_60000000"&
                                     genes_information$chr_pos < "5_61000000"),]$Gene)
 # Substitute
general_variants[which(general_variants$chr_pos == "5_60240986"),]$Gene = gene

# Search gene
gene = unique(genes_information[which(genes_information$chr_pos > "7_142400000"&
                                        genes_information$chr_pos < "7_142500000"),]$Gene)
# Substitute
general_variants[which(general_variants$chr_pos == "7_142457132"),]$Gene = gene


## OVARIAN GENE SEARCHING ##

sel = which(ovarian_variants$chr_pos %in% genes_information$chr_pos)# Select the common
for (i in sel){
  print(i)
  chrpos = ovarian_variants[i,]$chr_pos
  ovarian_variants[i,]$Gene = unique(genes_information[which(genes_information$chr_pos == chrpos),]$Gene)
}

# Check and manual fill
ovarian_variants[which(is.na(ovarian_variants$Gene)),]$chr_pos # "7_142457132"

# Search gene
gene = unique(genes_information[which(genes_information$chr_pos > "7_142400000"&
                                        genes_information$chr_pos < "7_142500000"),]$Gene)
# Substitute
ovarian_variants[which(ovarian_variants$chr_pos == "7_142457132"),]$Gene = gene


## BREAST GENE SEARCHING ##

sel = which(breast_variants$chr_pos %in% genes_information$chr_pos)# Select the common
for (i in sel){
  print(i)
  chrpos = breast_variants[i,]$chr_pos
  breast_variants[i,]$Gene = unique(genes_information[which(genes_information$chr_pos == chrpos),]$Gene)
}

# Check and manual fill
breast_variants[which(is.na(breast_variants$Gene)),]$chr_pos # "5_60240986"  "7_142457132"

# Search gene
gene = unique(genes_information[which(genes_information$chr_pos > "5_60200000"&
                                        genes_information$chr_pos < "5_60300000"),]$Gene)
# Substitute
breast_variants[which(breast_variants$chr_pos == "5_60240986"),]$Gene = gene

# Search gene
gene = unique(genes_information[which(genes_information$chr_pos > "7_142400000"&
                                        genes_information$chr_pos < "7_142500000"),]$Gene)
# Substitute
breast_variants[which(breast_variants$chr_pos == "7_142457132"),]$Gene = gene

## COLON GENE SEARCHING ##

sel = which(colon_variants$chr_pos %in% genes_information$chr_pos)# Select the common
for (i in sel){
  print(i)
  chrpos = colon_variants[i,]$chr_pos
  colon_variants[i,]$Gene = unique(genes_information[which(genes_information$chr_pos == chrpos),]$Gene)
}

# Check and manual fill
colon_variants[which(is.na(colon_variants$Gene)),]$chr_pos # "5_60240986"

# Search gene
gene = unique(genes_information[which(genes_information$chr_pos > "5_60200000"&
                                        genes_information$chr_pos < "5_60300000"),]$Gene)
# Substitute
colon_variants[which(colon_variants$chr_pos == "5_60240986"),]$Gene = gene


# install.packages('VennDiagram')
library(VennDiagram)

png("graphs/VennDiagram.png")
grid.newpage()
draw.triple.venn(area1 = length(unique(ovarian_variants$Gene)), # Ovarian
                 area2 = length(unique(breast_variants$Gene)),  # Breast
                 area3 = length(unique(colon_variants$Gene)), # Colon
                 n12 = length(intersect(ovarian_variants$Gene,breast_variants$Gene)), #Ovarian vs Breast
                 n23 = length(intersect(breast_variants$Gene,colon_variants$Gene)), #Breast vs Colon
                 n13 = length(intersect(ovarian_variants$Gene,colon_variants$Gene)),  # Ovarian vs Colon
                 n123 = length(intersect(intersect(ovarian_variants$Gene,breast_variants$Gene),colon_variants$Gene)), # ALL 
                 category = c("Ovarian", "Breast", "Colon"), # cancer types
                 lty = "dotted", lwd = 1, alpha = 0.3, 
                 cat.col = c("blue", "aquamarine3", "blueviolet"),
                 fill = c("blue", "aquamarine", "blueviolet"), scaled  = TRUE)
dev.off()

# Now we need to search each variant to determine the gene to which it belongs 

# Save the variants
general_variants %>% write.table("general_variants.txt",
                                 row.names = F,col.names = F,quote = F)
ovarian_variants %>% write.table("ovarian_variants.txt",
                                 row.names = F,col.names = F,quote = F)
breast_variants %>% write.table("breast_variants.txt",
                                row.names = F,col.names = F,quote = F)
colon_variants %>% write.table("colon_variants.txt",
                               row.names = F,col.names = F,quote = F)


#### TRANING AND TEST DATA ####

head(psam)

# Set seed to reproducibility
set.seed(123)

# GENERAL 
# Create the 70% random sample from psam (training set)
my_sample = sample(1:nrow(psam),round(nrow(psam)*0.70))

# Add the new columns 
psam$pheno_70 = NA
psam$pheno_30 = NA
# Fill the columns (number = sample used, NA = not used)
psam$pheno_70[my_sample] = psam$pheno[my_sample]
psam$pheno_30[-my_sample] = psam$pheno[-my_sample]

# OVARIAN
# Create the 70% random sample from psam (training set)
my_sample = sample(1:nrow(psam),round(nrow(psam)*0.70))

# Add the new columns 
psam$ovarian_70 = NA
psam$ovarian_30 = NA
# Fill the columns (number = sample used, NA = not used)
psam$ovarian_70[my_sample] = psam$ovarian[my_sample]
psam$ovarian_30[-my_sample] = psam$ovarian[-my_sample]

# BREAST
# Create the 70% random sample from psam (training set)
my_sample = sample(1:nrow(psam),round(nrow(psam)*0.70))

# Add the new columns 
psam$breast_70 = NA
psam$breast_30 = NA
# Fill the columns (number = sample used, NA = not used)
psam$breast_70[my_sample] = psam$breast[my_sample]
psam$breast_30[-my_sample] = psam$breast[-my_sample]

# COLON
# Create the 70% random sample from psam (training set)
my_sample = sample(1:nrow(psam),round(nrow(psam)*0.70))

# Add the new columns 
psam$colon_70 = NA
psam$colon_30 = NA
# Fill the columns (number = sample used, NA = not used)
psam$colon_70[my_sample] = psam$colon[my_sample]
psam$colon_30[-my_sample] = psam$colon[-my_sample]

# Write the psam 
write.table(psam,"../05_GWAS/outputs_plink/first_final_subset.psam",
            quote = F,row.names = F)
