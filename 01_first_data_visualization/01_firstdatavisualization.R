#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
###################### FIRST DATA VISUALIZATION ######################
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# We have the variants list from PANDORA

# We will evaluate this variants on the GCAT participants in order to calculate a PRS

# Libraries:

library(data.table) #fread
library(BiocGenerics) #distinct
library(dplyr)
library(xlsx)
library(openxlsx)
library(rgr)
library(ggplot2)
library(MASS)
library(plyr)
library(biomaRt)

setwd("/imppc/labs/dnalab/share/PRS/01_first_data_visualization/")

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
##################### Pandora SNPs vs SNP array ##################### 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #



# First we want to evaluate how many SNPs from the Pandora selection intersect with the ones in our SNP array

# Open the data we receive:

pandora_info <- read.xlsx("../Pandora_data/first_data/PandoraVariantsLibrary_20190129.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

head(pandora_info)

colnames(pandora_info)

# We open our curated SNP array:

SNParray_bim <- fread("../../data_storage/GCATCore_PLINK_curated/hg19/gcat_core.bim", head = F)

head(SNParray_bim)

length(unique(SNParray_bim$V2))

# We change the columns names

colnames(SNParray_bim) = c("chromosome", "id", "dist", "position", "ref", "alt")

head(SNParray_bim)

# Make a new column which contain the chr and pos in the same string

SNParray_bim$snps <- paste0("chr", SNParray_bim$chr, ":", SNParray_bim$pos)

# We want to obtain a similar file from pandora to be able to check both:

# delete the emptys rows
snps_pandora <- pandora_info$g..start 
snps_pandora = snps_pandora[-which(is.na(pandora_info$g..start))]
snps_pandora = as.data.frame(snps_pandora)


#to take the chromosome 
snps_pandora$chr = as.factor(gsub("chr", "", (tstrsplit(pandora_info$g..start[grep(":",pandora_info$g..start)],split=":")[[1]])))

#to take the position
snps_pandora$pos = as.integer(tstrsplit(pandora_info$g..start[grep(":",pandora_info$g..start)],split=":")[[2]])

head(snps_pandora)

colnames(snps_pandora) = c("snps", "chr", "pos")

length(unique(snps_pandora$snps))


# Now we make a left_join
SNParray_bim$chr = as.factor(SNParray_bim$chr)

union_snps = left_join(snps_pandora,SNParray_bim) %>% dplyr::select(snps,chr,pos,chromosome, position,id)

# The number of intersection snps
nrow(union_snps) - sum(is.na(union_snps$id))

head(union_snps)

sum(snps_pandora$snps %in% SNParray_bim$snps)

pandoravsGCAT = union_snps %>% filter(!is.na(id)) 

length(unique(pandoravsGCAT$id))

## In pandoravsGCAT we have the common snps for the array and pandora


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
################### Pandora SNPs vs imputed data #################### 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

#We have the pandora data 
head(snps_pandora)

# We need to take the imputed data:

imputed_data <- fread("../../data_storage/GCATCore_imputed_PLINK_august_18/info_all_genome.txt")

# See the distribution of the info of the intersection SNPs: SNP_id, MAF, info, reference_panel, chr_pos
head(imputed_data)
head(imputed_data$V5)

# the V5  is chr_position so lets add a new column in pandora with this info and set the same column name

snps_pandora$chr_pos = paste0(snps_pandora$chr, "_", snps_pandora$pos)

colnames(imputed_data) = c("V1", "MAF", "info", "ref_panel", "chr_pos")

dim(snps_pandora)
dim(imputed_data)
length(unique(imputed_data$chr_pos))


#Lets check the common SNPs

union_snps2 = left_join(snps_pandora,imputed_data) %>% dplyr::select(snps,chr,pos,chr_pos,MAF, info)

# The number of intersection snps
nrow(union_snps2) - sum(is.na(union_snps2$MAF))

#pandora vs imputed 
pandoravsimputed = union_snps2 %>% filter(!is.na(MAF)) 

length(unique(pandoravsimputed$chr_pos))

head(union_snps2)

png("graphs/PandoraSNPs_vs_GCATimputed.png")
hist(union_snps2$info, xlab = "GCAT imputed info", ylab = " ", main = "Pandora vs GCAT imputed")
dev.off()

# This imputed data have been filtered by info > 0.7 and maf of 0.01

# In pandora we have all the snps without filters

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
################### Pandora SNPs vs 1000 Genomes #################### 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# We want to extract from 1000Genomes the variants in order to make a PCA with the ICO-IDIBELL patients.

#Read the pvar file from 1000Genomes which contain the information of the variants  
milgen <- fread("../../data_storage/1000Genomes_Phase3_plink2/MAF_01/all_phase3.pvar", skip= 188,head = TRUE) #Which one
dim(milgen)
#We have Chromsome, position, id, ref, alt, qual, filter
head(milgen)

milgen$INFO = NULL
milgen$QUAL = NULL
milgen$FILTER = NULL


length(unique(milgen$ID)) # 14131756
# Make the intersect as in the previous step, to know which SNPs from pandora are present on milgenomes
#We have the pandora data 
head(snps_pandora)

#We will add a new column to the 1000G file to have chr_pos

milgen$chr_pos = paste0(milgen$`#CHROM`, "_", milgen$POS)
head(milgen)

#Now we make the union 

union_snps3 = left_join(snps_pandora,milgen) %>% dplyr::select(snps,chr,pos,chr_pos,ID)

head(union_snps3)

# The number of intersection snps
nrow(union_snps3) - sum(is.na(union_snps3$ID))

#pandora vs 1000G(MAF = 0.01) 
pandoravs1000G = union_snps3 %>% filter(!is.na(ID)) 

length(unique(pandoravs1000G$chr_pos))# 5182

#Filter the milgen data to get only the common variants 

# exclude different alleles in both files
length(union_snps3$ID)

common_snps = remove.na(union_snps3$ID, iftell = TRUE)[["x"]]

# Exclude duplicates 
length(common_snps)

common_snps = unique(common_snps)

common_snps %>% write.table("things/1000Gselected.txt",row.names = F,col.names = F,quote = F)


# With all this filters we update the data 
system(paste0("../../data_storage/1000Genomes_Phase3_plink2/plink2 ",
              "--pfile ../../data_storage/1000Genomes_Phase3_plink2/MAF_01/all_phase3 ",
              "--extract things/1000Gselected.txt --make-pgen --out things/1000Gselected")) 

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
############# 1000 Genomes + Pandora PCA   (MAF 0.01) ############### 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# we want to make a PCA with this SNPS to determinate the capacity of this SNPs to differenciate between populations

#Make the pca only with the SNPs selected from 1000Genomes
system(paste0("../../data_storage/1000Genomes_Phase3_plink2/plink2 ",
              "--pfile things/1000Gselected ",
              "--pca --threads 32 --out things/pandora_1000G_PCA"))

pvar = fread("things/1000Gselected.pvar", head = TRUE, skip = 188)

head(pvar)

psam = fread("things/1000Gselected.psam")

head(psam)


# read the eigenvec file 
pca_out = fread("things/pandora_1000G_PCA.eigenvec")

#We make a left join in order to have the population with the PCA data
pca_out = left_join(pca_out,psam)

head(pca_out)
str(pca_out)

# the population and superpopulation are factors
pca_out$Population = as.factor(pca_out$Population)
pca_out$SuperPop = as.factor(pca_out$SuperPop)

levels(pca_out$Population)
levels(pca_out$SuperPop)

#Set the colors as characters
pca_out$color = mapvalues(pca_out$SuperPop , from =levels(pca_out$SuperPop), 
                          to = c("darkmagenta", "darkred", "cyan3", "cornflowerblue", "darkgreen" ))
pca_out$color = as.character(pca_out$color)

head(pca_out)
table(pca_out$color)


#PLOT pcas
png("graphs/PandoraSNPs_1000Genomes.png",width = 16,height = 8,res = 300,units = "in")
par(mfrow=c(1,2))
plot(pca_out$PC1,pca_out$PC2,col=as.character(pca_out$color),asp=1,pch=16,cex=0.8,xlab="PC1",ylab="PC2")
legend("topright",col=c("darkmagenta", "darkred", "cyan3", "cornflowerblue", "darkgreen" ),
       c("AFR","AMR","EAS","EUR","SAS"), pch=16)
plot(pca_out$PC1,pca_out$PC3,col=as.character(pca_out$color),asp=1,pch=16,cex=0.8,xlab="PC1",ylab="PC3")
dev.off()

## Lets check the pca with all the 1000G snps

system(paste0("../../data_storage/1000Genomes_Phase3_plink2/plink2 ",
              "--pfile ../../data_storage/1000Genomes_Phase3_plink2/MAF_01/all_phase3 ",
              "--pca approx 3 --threads 32 --out things/1000G_PCA"))


# read the eigenvec file 
pca_all = fread("things/1000G_PCA.eigenvec")

#We make a left join in order to have the population with the PCA data
pca_all = left_join(pca_all,psam)

head(pca_all)
str(pca_all)

# the population and superpopulation are factors
pca_all$Population = as.factor(pca_all$Population)
pca_all$SuperPop = as.factor(pca_all$SuperPop)

levels(pca_all$Population)
levels(pca_all$SuperPop)

#Set the colors as characters
pca_all$color = mapvalues(pca_all$SuperPop , from =levels(pca_all$SuperPop), 
                          to = c("darkmagenta", "darkred", "cyan3", "cornflowerblue", "darkgreen" ))
pca_all$color = as.character(pca_all$color)

head(pca_all)
table(pca_all$color)


#PLOT pcas
png("graphs/1000G_PCA.png",width = 16,height = 8,res = 300,units = "in")
par(mfrow=c(1,2))
plot(pca_all$PC1,pca_all$PC2,col=as.character(pca_all$color),asp=1,pch=16,cex=0.8,xlab="PC1",ylab="PC2")
legend("bottomleft",col=c("darkmagenta", "darkred", "cyan3", "cornflowerblue", "darkgreen" ),
       c("AFR","AMR","EAS","EUR","SAS"), pch=16)
plot(pca_all$PC1,pca_all$PC3,col=as.character(pca_all$color),asp=1,pch=16,cex=0.8,xlab="PC1",ylab="PC3")
dev.off()

#We want to check the comparisson between pandora and all imputed,
#because in the previous case we use a imputation data fliter by maf 0.01

#We also want to know the population frequency of the pandora snps

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
################# 1000 Genomes + Pandora PCA all #################### 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

milgen_all <- fread("../../data_storage/1000Genomes_Phase3_plink2/all_phase3.pvar", skip= 188,head = TRUE) 

milgen_all$INFO = NULL 

head(milgen_all)

#Here we have the pandora SNPs
head(snps_pandora)
snps_pandora %>% write.table("data/snps_pandora.txt")
snps_pandora <- fread("data/snps_pandora.txt")
snps_pandora <- snps_pandora[,2:5]

#We will add a new column to the 1000G file to have chr_pos

milgen_all$chr_pos = paste0(milgen_all$`#CHROM`, "_", milgen_all$POS)
milgen_all %>% write.table("things/milgen_all.txt")
milgen_all <- fread("things/milgen_all.txt") 

dim(milgen_all)
head(milgen_all)

length(unique(milgen_all$ID))

#Now we make the union 

union_snps4 = left_join(snps_pandora,milgen_all) %>% dplyr::select(snps,chr,pos,chr_pos,ID)

union_snps4 %>% write.table("things/union_snps4.txt")
union_snps4 <- fread("things/union_snps4.txt")

tail(union_snps4)

#pandora vs 1000G(all) 
pandoravs1000Gall = union_snps4 %>% filter(!is.na(ID)) 

length(unique(pandoravs1000Gall$chr_pos))

# The number of intersection snps
nrow(union_snps4) - sum(is.na(union_snps4$ID))

#Filter the milgen data to get only the common variants 

# exclude different alleles in both files
length(union_snps3$ID)

common_snps = remove.na(union_snps3$ID, iftell = TRUE)[["x"]]

# Exclude duplicates 
length(snps_dup)
length(common_snps)

common_snps = unique(common_snps)

common_snps %>% write.table("things/1000Gselected_all.txt")
common_snps <- fread("things/1000Gselected_all.txt")

# With all this filters we update the data 
system(paste0("../../data_storage/1000Genomes_Phase3_plink2/plink2 ",
              "--pfile ../../data_storage/1000Genomes_Phase3_plink2/all_phase3 ",
              "--extract things/1000Gselected_all.txt --make-pgen --out things/1000Gselected_all")) 


#Make the pca only with the SNPs selected from 1000Genomes
system(paste0("../../data_storage/1000Genomes_Phase3_plink2/plink2 ",
              "--pfile things/1000Gselected_all ",
              "--pca --threads 32 --out things/pandora_1000G_PCA_all"))

pvar = fread("1000Gselected_all.pvar", head = TRUE, skip = 188)
dim(pvar)
head(pvar)
length(unique(pvar$ID))
psam = fread("things/1000Gselected_all.psam")

head(psam)


# read the eigenvec file 
pca_out = fread("things/pandora_1000G_PCA_all.eigenvec")

#We make a left join in order to have the population with the PCA data
pca_out = left_join(pca_out,psam)

head(pca_out)
str(pca_out)

# the population and superpopulation are factors
pca_out$Population = as.factor(pca_out$Population)
pca_out$SuperPop = as.factor(pca_out$SuperPop)

levels(pca_out$Population)
levels(pca_out$SuperPop)

#Set the colors as characters
pca_out$color = mapvalues(pca_out$SuperPop , from =levels(pca_out$SuperPop), 
                          to = c("darkmagenta", "darkred", "cyan3", "cornflowerblue", "darkgreen" ))
pca_out$color = as.character(pca_out$color)

head(pca_out)
table(pca_out$color)


#PLOT pcas
png("graphs/PandoraSNPs_1000Genomes.png",width = 16,height = 8,res = 300,units = "in")
par(mfrow=c(1,2))
plot(pca_out$PC1,pca_out$PC2,col=as.character(pca_out$color),asp=1,pch=16,cex=0.8,xlab="PC1",ylab="PC2")
legend("bottomright",col=c("darkmagenta", "darkred", "cyan3", "cornflowerblue", "darkgreen" ),
       c("AFR","AMR","EAS","EUR","SAS"), pch=16)
plot(pca_out$PC1,pca_out$PC3,col=as.character(pca_out$color),asp=1,pch=16,cex=0.8,xlab="PC1",ylab="PC3")
dev.off()

dim(pca_out)

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
################## pop freq of the common snps ###################### 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# We download Output snp151 as a BED table in order to check the pandora SNPs on it
# NCBI dbSNP Build 151 (Oct 06, 2017)
#Number of (ss#'s)  with frequency 364,060,923
#Number of (rs#'s) in gene 381,785,470
#Number of Submissions (ss#'s) 1,803,563,957

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
####################### dbSNP vs Pandora ########################### 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# We download Output snp151 as a BED table in order to check the pandora SNPs on it
# NCBI dbSNP Build 151 (Oct 06, 2017)
#Number of (ss#'s)  with frequency 364,060,923
#Number of (rs#'s) in gene 381,785,470
#Number of Submissions (ss#'s) 1,803,563,957

# Read the dbSNP table 
dbSNP_variation<- fread("../../imppc/labs/dnalab/share/dbSNP_variation", head = F, fill=TRUE)

#Repeat 4 times
dbSNP_variation = dbSNP_variation[-nrow(dbSNP_variation),]
dbSNP_variation = dbSNP_variation[-nrow(dbSNP_variation),]
dbSNP_variation = dbSNP_variation[-nrow(dbSNP_variation),]
dbSNP_variation = dbSNP_variation[-nrow(dbSNP_variation),]

head(dbSNP_variation)
tail(dbSNP_variation)
dim(dbSNP_variation)

# Change the col names 
colnames(dbSNP_variation) <- c("chromosome", "start", "end", "id", "V5", "strand")
head(dbSNP_variation)

# Delete the chr from the chromosome column
dbSNP_variation$chromosome = as.factor(gsub("chr", "", dbSNP_variation$chromosome))
head(dbSNP_variation)

# Create the chromosome + position column
dbSNP_variation$chr_pos = paste0(dbSNP_variation$chromosome, "_", dbSNP_variation$start)
head(dbSNP_variation)

#Now read the pandora's snps
snps_pandora <- fread("data/snps_pandora.txt")
snps_pandora <- snps_pandora[,2:5] #Delete the first column
head(snps_pandora)
tail(snps_pandora)

# Now we have both with a chr_pos column lets make a left join
union_snps5 = left_join(snps_pandora,dbSNP_variation) %>% dplyr::select(snps,chr,pos,id)
head(union_snps5)

# The number of intersection snps:
nrow(union_snps5) - sum(is.na(union_snps5$id))

dim(dbSNP_variation)
dim(snps_pandora)

#Pandora vs dbSNP
pandoravsdbSNP = union_snps5 %>% filter(!is.na(id)) 
length(unique(pandoravsdbSNP$id))

# Histogram of the freq inside pandora
# we make this transformation because if there is a freq higher than 0.5 its not the maf
pandora_info$inHouseFrequency = ifelse(pandora_info$inHouseFrequency>0.5,
                                       1-pandora_info$inHouseFrequency,
                                       pandora_info$inHouseFrequency)

# Histogram of the frequencies 
png("graphs/Inhousefreq_histogram.png")
hist(pandora_info$inHouseFrequency, main = "histogram Pandora", 
     xlab = "Pandora in house frequency", ylab = "", breaks = 40)
abline(v = 0.001, col = "darkgreen", lwd=5, lty=3)
abline(v = 0.01, col = "darkred", lwd=5, lty=3)
abline(v = 0.1, col = "cyan3", lwd=5, lty=3)
legend("topright",col=c("darkgreen", "darkred", "cyan3"),
       c("maf = 0.001","maf = 0.01","maf = 0.1"), pch=16)
dev.off()

dim(pca_out)

# Check the number of snps under this frequencies 
length(which(pandora_info$inHouseFrequency<0.1))
length(which(pandora_info$inHouseFrequency<0.01))
length(which(pandora_info$inHouseFrequency<0.001))
length((snps_pandora$chr_pos))

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
################### Pandora SNPs in HRC panel ####################### 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# First we obtain a list of the files (with the path) inside the HRC raw data 
my_files = list.files("../../data_storage/GCATCore_imputed_raw_august_18/GCAT imputation August 2018/HRC/",recursive = T)

# We filter all the files of this list because we are only interested in the info files
my_files = my_files[grep(".impute_info",my_files)]
length(my_files)
head(my_files)

#We create a directory on the wdir called HRC
dir.create("HRC")

# Make a loop to search all this files and look for the pandora's snps

for(i in 1:length(my_files)){
  # This line is in order to avoid errors when it find an empty file 
  if(length(readLines(paste0("../../data_storage/GCATCore_imputed_raw_august_18/GCAT imputation August 2018/HRC/",
                             my_files[i])))!=0){
    # info_file is the path of the current file 
    info_file = fread(paste0("../../data_storage/GCATCore_imputed_raw_august_18/GCAT imputation August 2018/HRC/",
                             my_files[i]))
    # obtain the chromosome 
    chr = tstrsplit(my_files[i],"/")
    chr = gsub("Chr_","",chr[[2]])
    # create the chromosome_position
    info_file$chr_pos = do.call(paste0,list(chr,"_",info_file$position))
    # If there is some common snp
    if(sum(info_file$chr_pos %in% snps_pandora$chr_pos)!=0){
      # Obtain the ones with correspondence between both 
      info_file = info_file[info_file$chr_pos %in% snps_pandora$chr_pos,] 
      # Write the common snps in a txt file (apped in each iteration of the loop)
      fwrite(info_file,"HRC/info_file.txt",row.names = F,col.names = F,append = T)
    }
    #This print show the progress  
    print(i)
  }
}

# Read the info file 
info_file = fread("HRC/info_file.txt")
head(info_file)

# SNPs in pandora and raw data = 6779
dim(info_file)

hist(info_file$V6, main = "HRC and pandora common snps", ylab = "", xlab = "MAF") # histogram of maf
hist(info_file$V7, main = "HRC and pandora common snps", ylab = "", xlab = "info") # histogram of the info 

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
################## Pandora SNPs in all panels ####################### 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# Find the files in all the panels and then repeat the same as in HRC
my_files1 = list.files("../../data_storage/GCATCore_imputed_raw_august_18/GCAT imputation August 2018/1kgphase3/",recursive = T)
my_files2 = list.files("../../data_storage/GCATCore_imputed_raw_august_18/GCAT imputation August 2018/gonl/",recursive = T)
my_files3 = list.files("../../data_storage/GCATCore_imputed_raw_august_18/GCAT imputation August 2018/HRC/",recursive = T)
my_files4 = list.files("../../data_storage/GCATCore_imputed_raw_august_18/GCAT imputation August 2018/uk10k/",recursive = T)

length(my_files1)
length(my_files2)
length(my_files3)
length(my_files4)

my_files1 = my_files1[grep(".impute_info",my_files1)]
my_files2 = my_files2[grep(".impute_info",my_files2)]
my_files3 = my_files3[grep(".impute_info",my_files3)]
my_files4 = my_files4[grep(".impute_info",my_files4)]

my_files1 = paste0("1kgphase3/", my_files1)
my_files2 = paste0("gonl/", my_files2)
my_files3 = paste0("HRC/", my_files3)
my_files4 = paste0("uk10k/", my_files4)

length(my_files1)
length(my_files2)
length(my_files3)
length(my_files4)

my_files = c(my_files1, my_files2, my_files3, my_files4)
length(my_files)

dir.create("all_panels")

for(i in 1:length(my_files)){
  
  if(length(readLines(paste0("../../data_storage/GCATCore_imputed_raw_august_18/GCAT imputation August 2018/",
                             my_files[i])))!=0){
    
    info_file = fread(paste0("../../data_storage/GCATCore_imputed_raw_august_18/GCAT imputation August 2018/",
                             my_files[i]))
    
    chr = tstrsplit(my_files[i],"/")
    panel = chr[[1]]
    chr = gsub("Chr_","",chr[[3]])
    
    #Add a new variable (col) named panel to indicate the procedence of the information
    info_file$panel = panel
    
    info_file$chr_pos = do.call(paste0,list(chr,"_",info_file$position))
    
    if(sum(info_file$chr_pos %in% snps_pandora$chr_pos)!=0){
      
      info_file = info_file[info_file$chr_pos %in% snps_pandora$chr_pos,] 
      
      fwrite(info_file,"all_panels/info_file_all.txt",row.names = F,col.names = F,append = T)
    }
    
    print(i)
  }
}


# SNPs in pandora and raw data 12772
info_file2 = fread("all_panels/info_file_all.txt")
colnames(info_file2) <- colnames(info_file)
head(info_file2)

# We see the dimension but here we have information from the 4 panels and there are repetitions
dim(info_file2)

# See the unique snps 
unique_snps <- unique(info_file2$chr_pos)
length(unique_snps)

# Now we want to obtain the best panel for each snp, the one with higher info
df = as.data.table(info_file2)
df_intersection = df[df[, .I[which.max(info)], by = chr_pos]$V1,]

#Save the results:
fwrite(df_intersection,"all_panels/intersection_snps.txt",row.names = F,col.names = T)

#Read the data 
df_intersection = fread("all_panels/intersection_snps.txt")
dim(df_intersection)
head(df_intersection)

#Number of SNPs for each panel 
table(df_intersection$panel)

#Histogram of INFO
png("graphs/info_4panels.png")
hist(df_intersection$info, main = "Info Pandora vs 4th panels", 
     xlab = "SNPs info", ylab = "", breaks = 20)
abline(v = 0.7, col = "darkgreen", lwd=5)
abline(v = 0.1, col = "darkred", lwd=5)
legend("top",col=c("darkgreen", "darkred"),
       c("info = 0.7","info = 0.1"), pch=16)
dev.off()

sum(df_intersection$info > 0.7)
sum(df_intersection$info < 0.1)

#Histogram of the maf

df_intersection$exp_freq_a1 = ifelse(df_intersection$exp_freq_a1>0.5,
                                     1-df_intersection$exp_freq_a1,
                                     df_intersection$exp_freq_a1)

png("graphs/maf_4panels.png")
hist(df_intersection$exp_freq_a1, xlab= "maf"
     , ylab = "", main = "Pandora SNPs commmon with the 4 panels", breaks = 50)
abline(v = 0.001, col = "darkred", lwd=5)
legend("top",col=c("darkred"),
       c("maf = 0.001"), pch=16)
dev.off()

sum(df_intersection$exp_freq_a1 == 0)

hist(df_intersection[which(df_intersection$exp_freq_a1 == 0)]$info, 
     main = "Monomorphic common site", xlab = "info", ylab = "" )

mono <- df_intersection[which(df_intersection$exp_freq_a1 == 0)]

mono_high_info <- mono[which((mono$info) > 0.7)]
head(mono_high_info)
dim(mono_high_info)

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
######################## Check Pandora ############################## 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# Check how many snps with rs we have in pandora = 15546
length(pandora_info$dbSNP)
pan = pandora_info %>% filter(!is.na(dbSNP)) %>% dplyr::select(g..start,Ref,Alt,inROI, dbSNP)
head(pan)
dim(pan)

head(df_intersection)

# we need to obtain chr_pos in pan

#to take the chromosome 
pan$chr = as.factor(gsub("chr", "", (tstrsplit(pan$g..start[grep(":",pan$g..start)],split=":")[[1]])))

#to take the position
pan$pos = as.integer(tstrsplit(pan$g..start[grep(":",pan$g..start)],split=":")[[2]])

pan$chr_pos = do.call(paste0,list(pan$chr, "_", pan$pos))

head(pan)


# Left join
union_snps6 = left_join(df_intersection,pan) %>% dplyr::select(chr,pos,chr_pos,inROI, dbSNP, Ref, Alt)
head(union_snps6)


# The number of intersection snps
nrow(union_snps6) - sum(is.na(union_snps6$dbSNP))
tail(union_snps6)
dim(df_intersection)
dim(pan)

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
###################### Genes in Pandora ############################# 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# genes in pandora
table(pandora_info$Gene)

# Result: 140 genes 
pandora_genes <- unique(pandora_info$Gene[!is.na((pandora_info$Gene))])
length(pandora_genes)

png("graphs/genes_in_pandora.png")
hist(table(pandora_info$Gene), main = "Pandora genes (140)", ylab= "", xlab = "number of SNPs per gene")
dev.off()


#lets clean pandora info
pandora_info <- read.xlsx("../Pandora_data/first_data/PandoraVariantsLibrary_20190129.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
pandora_info = pandora_info %>% dplyr::select(Gene, transcriptRegion, g..start, Ref, Alt, dbSNP, inHouseFrequency, inHouseNumber, popFreqMax, n_1000g_eur)
head(pandora_info)
head(snps_pandora)

dim(pandora_info)
dim(snps_pandora)

tail(pandora_info)
tail(snps_pandora)

genes = as.data.frame(pandora_genes)
chr = c()
min_pos = c()
max_pos = c()
g = "EXO1"
for (g in pandora_genes){
  #Select the rows of the gene
  sel = pandora_info %>% filter(Gene == g) 
  #For chromosome
  chromosome = gsub("chr", "", tstrsplit(sel[1,]$g..start,split = ":")[[1]])
  chr <- c(chr, chromosome)
  #For position
  mini = tstrsplit(sel$g..start,split = ":")[[2]] %>% min()
  maxi = tstrsplit(sel$g..start,split = ":")[[2]] %>% max()
  min_pos <- c(min_pos, mini)
  max_pos <- c(max_pos, maxi)
}

#We have here the genes and che chromosome for each gene
genes$chr = chr
genes$min_snp = min_pos
genes$max_snp = max_pos### The min and max position are the min or max snp position find


head(genes)

## Know we want to study the distance of this genes 
#for that we need to search the higher and lower position of the rois
#for each gene filter the chromosome in rois then the position and make max or min
#g = "A2ML1"
min_pos = c()
max_pos = c()
g = "EXO1"
for (g in pandora_genes){
  #Take the chromosome
  chr = (genes[genes$pandora_genes == g,])$chr
  min = as.integer((genes[genes$pandora_genes == g,])$min_snp)
  max = as.integer((genes[genes$pandora_genes == g,])$max_snp)
  #filter the rois
  # all the rois in this chromosome
  selection <- (rois[rois$V1 == chr,])
  #for the start:
  #we need to filter all the starts smaller than the first SNP
  start <- selection %>% filter(selection$V2 <= min)
  #Then select the maximun 
  start <- max(as.vector(start$V2))
  #for the end:
  #we need to filter all the ends higher than the last SNP
  end <- selection %>% filter(selection$V3 >= max)
  #Then select the minimun
  end <- min(as.vector(end$V3))
  min_pos <- c(min_pos, start)
  max_pos <- c(max_pos, end)
}

length(min_pos)
length(max_pos)

genes$start = min_pos
genes$end = max_pos

head(genes)
dim(genes)

#PROBLEM some snps start is smaller than the smaller
#start roi and the same for end 
#Change the -inf for the minsnp
#Change the inf for the maxsnp
gene = "FANCD2"
for (gene in pandora_genes){
  if (genes[genes$pandora_genes == gene,]$start == "-Inf"){
    genes[genes$pandora_genes == gene,]$start = genes[genes$pandora_genes == gene,]$min_snp
  }
  if (genes[genes$pandora_genes == gene,]$end == "Inf"){
    genes[genes$pandora_genes == gene,]$end = genes[genes$pandora_genes == gene,]$max_snp
  }  
}


pandora_info[47343,]

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
####################### Groups of Genes ############################# 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

head(genes)

colnames(genes) = c("name", "chr", "first snp", "last snp", "start", "end")

head(genes)


# check the overlap of the genes
genesforoverlap = genes
overlap = c()
gene = 1
comp = 1
comp = 138
for (gene in 1:139){#from 1 to 140 1 by 1
  #For each of the genes we define the start and end
  start <- genes[gene]$start
  end <- genes[gene]$end
  #Now we need to check the other genes 
  for (comp in 2:nrow(genesforoverlap)){
    # start position in interval
    if (genesforoverlap[comp]$start > start & genesforoverlap[comp]$start < end){
      towrite = genesforoverlap[comp]
      towrite$overlap = genesforoverlap[gene]$name
      overlap = rbind(overlap, towrite)
      break
    }
    if (genesforoverlap[comp]$end > start & genesforoverlap[comp]$end < end){
      towrite = genesforoverlap[comp]
      towrite$overlap = genesforoverlap[gene]$name
      overlap = rbind(overlap, towrite)
    }
  }
  #Delete the first gene
  genesforoverlap = genesforoverlap[2:nrow(genesforoverlap),]
}

dim(overlap)
head(overlap)


#Let's plot the intervals in order to see the overlap
#One plot for each chromosome

all_chr = c(1:19)
all_chr = as.character(all_chr)
all_chr = c(all_chr,c("22", "X"))

chromosome = 4
chromosome = 1

for (chromosome in all_chr){
  g = genes[genes$chr == chromosome,]
  dat <- data.frame(pos = c(1:nrow(g)),
                    start = as.vector(g$start),
                    end = as.vector(g$end),
                    names = as.vector(g$name))
  plot(c(min(g$start),max(g$end)),c(1,(nrow(g)+0.2)),type = "n",axes = FALSE,xlab = "Position",
       ylab = "",main=(paste0("Chromosome ", chromosome, ": ROIs positions.")))
  segments(x0 = dat$start,
           y0 = dat$pos,
           x1 = dat$end,
           y1 = dat$pos,
           col = "blue",
           lwd = 6,
           lend = 10)
  text(x = ((dat$start)+10000000),y = ((dat$pos)+ 0.15), labels = dat$names,font = 0.5)
  axis(1,at = c((as.integer(min(g$start))),(as.integer(max(g$end)))),
       labels = TRUE,tcl = 0.5, pos = -1)
}

#Lets try to make the same plot but with start snps and end snp CHANGE

for (chromosome in all_chr){
  g = genes[genes$chr == chromosome,]
  dat <- data.frame(pos = c(1:nrow(g)),
                    start = as.vector(g$`first snp`),
                    end = as.vector(g$`last snp`),
                    names = as.vector(g$name))
  plot(c(min((g$`first snp`)+0.2),max((g$`last snp`)+0.2)),c(1,(nrow(g))+0.2),type = "n",axes = FALSE,
       xlab = "Position",ylab = "",main=(paste0("Chromosome ", chromosome, ": SNPs positions.")))
  segments(x0 = dat$start,
           y0 = dat$pos,
           x1 = dat$end,
           y1 = dat$pos,
           col = "blue",
           lwd = 6,
           lend = 10)
  text(x = ((dat$start)+10000000),y = ((dat$pos)+ 0.15), labels = dat$names,font = 0.5)
  axis(1,at = c((as.integer(min(g$start))),(as.integer(max(g$end)))),labels = TRUE,tcl = 0.5, pos = -1)
}

#Problem.  -inf inf some snps out of rois boundaries

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
################### Groups of Genes in ROIs ######################### 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

#Now we need in ROI
pandora_info <- read.xlsx("../Pandora_data/first_data/PandoraVariantsLibrary_20190129.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
pandora_info = pandora_info %>% dplyr::select(Gene, transcriptRegion, g..start, Ref, Alt, inROI, dbSNP, inHouseFrequency, inHouseNumber, popFreqMax, n_1000g_eur)
head(pandora_info)

min_pos_inROI = c()
max_pos_inROI = c()
g = "EXO1"
for (g in pandora_genes){
  #Select the rows of the gene with the IN ROI filter
  sel = pandora_info %>% filter(Gene == g) %>% filter(inROI == TRUE)
  #For chromosome
  chromosome = gsub("chr", "", tstrsplit(sel[1,]$g..start,split = ":")[[1]])
  #For position
  mini = tstrsplit(sel$g..start,split = ":")[[2]] %>% min()
  maxi = tstrsplit(sel$g..start,split = ":")[[2]] %>% max()
  min_pos_inROI <- c(min_pos_inROI, mini)
  max_pos_inROI <- c(max_pos_inROI, maxi)
}

head(genes)

#Add in first snp and last snps in ROI variables
genes$`first snp ROI` = min_pos_inROI
genes$`last snp ROI` = max_pos_inROI

head(genes)

#Now we want to obtain as in the previous case the start and end
#But taking into acc the roi

min_pos = c()
max_pos = c()
g = "PPP1CB"
for (g in pandora_genes){
  #Take the chromosome
  chr = (genes[genes$name == g,])$chr
  min = as.integer((genes[genes$name == g,])$`first snp ROI`)
  max = as.integer((genes[genes$name == g,])$`last snp ROI`)
  #filter the rois
  # all the rois in this chromosome
  selection <- (rois[rois$V1 == chr,])
  #for the start:
  #we need to filter all the starts smaller than the first SNP
  start <- selection %>% filter(selection$V2 <= min)
  #Then select the maximun 
  start <- max(as.vector(start$V2))
  #for the end:
  #we need to filter all the ends higher than the last SNP
  end <- selection %>% filter(selection$V3 >= max)
  #Then select the minimun
  end <- min(as.vector(end$V3))
  min_pos <- c(min_pos, start)
  max_pos <- c(max_pos, end)
}


################ TAKE CARE!
# There are min snps inside inROI which are smaller than the smaller roi start in this chr

#> min (this is the smaller snp in roi inside the PPP1CB gene)
#[1] 28975017
#> min(selection$V2) #from all the rois in this chromosome the smaller
#[1] 29416069


#32890563   

length(min_pos)
length(max_pos)

genes$start = min_pos
genes$end = max_pos



#When we select for example for gene EXO1 without in roi selection we have 408, with in roi selection 107

#Change the chromosome to integer to be able to sort 
genes$chr = as.integer(genes$chr)
str(genes)
#Sort by chromosome
genes <- genes[order(genes$chr),]

#Sustitute the NAs by X
genes[is.na(genes)] <- "X"

#save the genes
fwrite(genes,"data/pandora_genes.txt")

genes <- fread("data/pandora_genes.txt")


#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
################# pop freq of the common snps ####################### 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

head(pandora_info)

# Why we only find 1/3 of the snps in all 1000genomes data?

pandora_info[,pandora_info$g..start == "chr1:45794801"]

rois <- fread("../Pandora_data/first_data/rois/rois.bed")

head(rois)
str(rois)
colnames(rois) <- c("chr", "start", "end", "V4", "V5", "V6")

rois %>% filter(chr == "1") %>% dplyr::select(start) %>% min()
rois %>% filter(chr == "1") %>% dplyr::select(end) %>% max()

str(pandora_info)
length(which(pandora_info$inROI)) #TRUE
length(which(!pandora_info$inROI)) #FALSE

length(pandora_info$inROI)

head(snps_pandora)
head(pandora_info)

head(df_intersection)
dim(snps_pandora) #unique = 45489
dim(df_intersection) #unique = 17252

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
############## Check the SNPs not found in panels ################### 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

# The SNPs we DONT have in the panels
#First step check the snps not in the panels
notin_panel <- snps_pandora[!(snps_pandora$chr_pos %in% df_intersection$chr_pos),]
head(notin_panel)
dim(notin_panel)
#take the info from pandora_info for all this not in panel snps
notin_panel <- pandora_info[pandora_info$g..start %in% notin_panel$snps,]
head(notin_panel)
dim(notin_panel)

#Check the length of the unique not in panel snps and filter to get only unique
length(unique(notin_panel$g..start))#28237
dim(notin_panel)
head(notin_panel)

length(which(notin_panel$inHouseFrequency<0.1))
length(which(notin_panel$inHouseFrequency<0.01))
length(which(notin_panel$inHouseFrequency<0.001))

sum(is.na(pandora_info$g..start))

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
#################### Repeated pandora SNPs ########################## 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

head(pandora_info)

#This are the unique pandora_snps by chr:position
unique_snps <- unique(pandora_info$g..start)

#total pandora length is 47569
dim(pandora_info)

#The total number of pandora SNPs without g..start is 197
sum(is.na(pandora_info$g..start))

unique_snps <- sort(unique_snps)
head(unique_snps)
tail(unique_snps)

duplicated_snps = c()

snp= "chr12:112893676"
for (snp in unique_snps){
  sel <- pandora_info[which(pandora_info$g..start == snp),]
  if (nrow(sel) > 1){
    duplicated_snps = rbind(duplicated_snps, sel)
    print(snp)
  }
}

str(duplicated_snps)

#duplicated_snps <- duplicated_snps %>% dplyr::select(Gene, NM, transcriptRegion, g..start, Ref, Alt, 
#                      inROI, reasoning, dbSNP, clinvar, inHouseFrequency,
#                      inHouseNumber, popFreqMax, n_1000g_eur, exac_nfe)

fwrite(duplicated_snps,"things/pandora_duplicated.txt")

duplicated_snps2 <- duplicated_snps %>% dplyr::select(Gene, g..start, Ref, Alt)

#Now we have all the duplicated positions 
#but a fast check shows that these are different insertions, deletions duplications or sustitutions
#lets check if they are exactly the same

real_dup = c()
snp= "chr11:108151707"
for (snp in unique(duplicated_snps2$g..start)){
  sel <- duplicated_snps2[which(duplicated_snps2$g..start == snp),]
  sel2 <- sel[which(duplicated(sel)),]
  if ((nrow(sel2)) > 0){
    print(snp)
    sel3 <- duplicated_snps %>% filter(duplicated_snps$g..start == snp)
    sel3 <- sel3 %>% filter(sel3$Alt == sel2$Alt)
    sel3 <- sel3 %>% filter(sel3$Ref == sel2$Ref)
    real_dup = rbind(real_dup, sel3)
  }
}

head(real_dup)
fwrite(real_dup,"things/pandora_REAL_duplicated.txt")

real_dup <- fread("things/pandora_REAL_duplicated.txt")

real_dup[19:20,] # SHOW


write.csv(real_dup, file = "things/dup_pandora_snps.csv")

#Now we want to obtain all the info from that snps 

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
######################## 800 risk SNPs ############################## 
#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

#OjO in order to read this we need to use the NOT console Rstudio version
riskSNPs <- readWorksheetFromFile("data/823_SNPs_TCGA-type.xlsx", 
                                  sheet=2, 
                                  startRow = 1,
                                  endCol = 9)
head(riskSNPs)
head(pandora_info$g..start)

fwrite(riskSNPs,"riskSNPs.txt")
riskSNPs <- fread("riskSNPs.txt")

# In risk we have the chr and position of each snp 
riskSNPs$g..start = paste0("chr",riskSNPs$chr_name, ":", riskSNPs$chrom_start)
head(riskSNPs)

# Check pandora info
str(pandora_info)

#Select only the interesting information
risk_pandora <- pandora_info %>% dplyr::select(Gene, NM, transcriptRegion, g..start, Ref, Alt, 
                                               inROI, reasoning, dbSNP, clinvar, inHouseFrequency,
                                               inHouseNumber, popFreqMax, n_1000g_eur, exac_nfe)
head(risk_pandora)
head(riskSNPs)

#Now we make a left-join to see the correspondence

risk_pandora <- left_join(risk_pandora, riskSNPs) %>% dplyr::select(Gene, NM, transcriptRegion, Chromosome, chrom_start, g..start, Ref, Alt, inROI, dbSNP,
                                                                    clinvar, inHouseFrequency, inHouseNumber, popFreqMax, n_1000g_eur, exac_nfe,
                                                                    SNP, Source, TCGA.cancer, allele)
head(risk_pandora)

#Now letf check

# As we can see there are 5 SNPs in common between the 800 of risk and pandora

risk_pandora_snps <- risk_pandora[!is.na(risk_pandora$SNP),]
nrow(risk_pandora_snps)

