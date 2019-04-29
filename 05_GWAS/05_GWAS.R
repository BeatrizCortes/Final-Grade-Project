######################################################################################
# 05 GWAS Preparation                                                                #
######################################################################################

# Libraries
library(openxlsx)
library(dplyr)
library(data.table)

# Set the working directory 
setwd("/imppc/labs/dnalab/share/PRS/05_GWAS/")
#dir.create("outputs_plink")

# With PLINK we obtain the frequence of the prs_data: all the data Pandora+GCAT together
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.psam ",
              " --freq ",
              " --out /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/prs_data"))

# Read the frequences 
freq = fread("outputs_plink/prs_data.afreq")

head(freq)

# Create the MAF variable in freq
freq$MAF = ifelse(freq$ALT_FREQS>0.5,1-freq$ALT_FREQS,freq$ALT_FREQS)

# Check the maf with a histogram
png("graphs/MAF_histogram.png")
hist(freq$MAF, main = "", 
     xlab = "MAF", ylab = "Frequency", breaks = 23)
dev.off()

# As we can see in general we have small MAF, rare variants

# Now we make a PCA with a filter to the rare variants (maf > 0.01)
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.psam ",
              " --pca approx --maf 0.01 ",
              " --out /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/prs_data"))

# Read the Pca eigenvectors values 
pca_out = fread("outputs_plink/prs_data.eigenvec")

head(pca_out)
dim(pca_out)

# Prepare the colours
my_cols = c(rep(1,4988),rep(2,1581))

# Make the plot
png("graphs/pca.png")
plot(pca_out$PC1,pca_out$PC2,pch=16,col=my_cols, 
     main = "Pandora vs GCAT PCA", 
     xlab= "PC1", ylab = "PC2")
legend("bottomleft",legend = c("GCAT", "Pandora"), col = c("black", "red"), 
       title = "Populations", pch = 16)
dev.off()

# The result is horrible, as we can see the two groups (Pandora and GCAT) are totally separated 

# GCAT-Pandora freqs #

# Read the data
my_sample = fread("prs_data.psam")

# Divide the data into the gcat and the pandora individuals
my_sample[1:4988,] %>% write.table("outputs_plink/gcat.txt",row.names = F,quote = F)
my_sample[4989:6569,] %>% write.table("outputs_plink/pandora.txt",row.names = F,quote = F)

# GCAT: now calculate the frequencies only for GCAT

system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.psam ",
              " --freq ",
              " --keep outputs_plink/gcat.txt",
              " --out /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/gcat"))

# Read the frequencies of GCAT
freq_g = fread("outputs_plink/gcat.afreq")

# Change the names to have different names and be able to make the merge with pandora data
colnames(freq_g)[5:6] = c("ALT_FREQS_GCAT","OBS_CT_GCAT")


# Pandora: repeat the same process with the pandora data 
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.psam ",
              " --freq ",
              " --keep outputs_plink/pandora.txt",
              " --out /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/pandora"))

freq_p = fread("outputs_plink/pandora.afreq")
colnames(freq_p)[5:6] = c("ALT_FREQS_PANDORA","OBS_CT_PANDORA")

# Now we can make the merge of the data 
freqs = left_join(freq_g,freq_p)

head(freqs)
colnames(freqs) = c("chr", 
                    "ID", 
                    "Ref_freqs", 
                    "Alt_freqs",
                    "ALT_FREQS_GCAT", 
                    "N_GCAT", 
                    "ALT_FREQS_PANDORA", 
                    "N_Pandora")
head(freqs)

# Now we calculate the difference for each frequency between Pandora and GCAT
freqs$diff_Pandora_GCAT = abs(freqs$ALT_FREQS_GCAT-freqs$ALT_FREQS_PANDORA)

# We order the variants in function of this difference to check the higher differences
freqs %>% arrange(desc(diff_Pandora_GCAT)) %>% head()

# As we can see the differences are high

# Read the Pandora information
variants_library = read.xlsx("../02_Pandora_data_preparation/first_data/PandoraVariantsLibrary_20190129.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
colnames(variants_library)[11] = "gstart" #Change the name of this column
dim(variants_library)

# Delete the empty gstart variants
variants_library = variants_library[!is.na(variants_library$gstart),]

# To take the chromosome 
variants_library$chr = as.factor(gsub("chr", "", (tstrsplit(variants_library$gstart[grep(":",variants_library$gstart)],split=":")[[1]])))

# To take the position
variants_library$pos = as.integer(tstrsplit(variants_library$gstart[grep(":",variants_library$gstart)],split=":")[[2]])

# To create the chr_pos
variants_library$chr_pos = paste0(variants_library$chr,"_",variants_library$pos)

# Read the bim file of prs_data
bim = fread("../05_GWAS/prs_data.bim")

# Create also the chr_pos
bim$chr_pos = paste0(bim$V1,"_",bim$V4)

# add to freqs chr_pos from bim
freqs$chr_pos = bim$chr_pos

# create the pandora_freqs taking only the variants of variants_library present in the bim file
pandora_freqs = variants_library[variants_library$chr_pos %in% bim$chr_pos,]

# Now add the information of the variants_library to the freqs 
freqs = left_join(freqs,(variants_library %>% select(chr_pos,Ref,Alt,inHouseNumber,inHouseFrequency,n_1000g_eur)))

# We have: 
# Chrosome, ID, REF, ALT, Alt_freq_GCAT, Observations_GCAT, Alt_freqs_Pandora, Observations_Pandora, diff_pandora_GCAT and chr_pos, Ref, Alt, inHourse_Number, inHouseFrequency, N_1000g_eur
head(freqs)

# Create a subset only with the variants with correspondent alleles in both datasets
freqs_subset = freqs[which((freqs$REF==freqs$Ref & freqs$ALT==freqs$Alt) | 
            (freqs$REF==freqs$Alt & freqs$ALT==freqs$Ref)), ]

# Now we check the changed 
which_different = which(freqs_subset$REF!=freqs_subset$Ref)

# And change the frequencies
freqs_subset$n_1000g_eur[which_different] = 1-freqs_subset$n_1000g_eur[which_different] 

head(freqs_subset)

# Now we make the difference between GCAT and 1000G and we see the differences are very small
# This make us see that this big differences between GCAT and Pandora are not a GCAT problem
freqs_subset$diff_gcat_1000G = abs(freqs_subset$ALT_FREQS_GCAT-freqs_subset$n_1000g_eur)

png("graphs/differences_GCAT_1000G.png")
hist(freqs_subset$diff_gcat_1000G, main = "GCAT vs 1000Genomes",
     xlab = "Frequency differences", ylab = "")
dev.off()

png("graphs/differences_GCAT_Pandora.png")
hist(freqs_subset$diff_Pandora_GCAT, main = "GCAT vs Pandora",
     xlab = "Frequency differences", ylab = "")
dev.off()

# Make the subset and the plot
hist_subset = freqs_subset %>% filter(diff_Pandora_GCAT >0.1)
png("graphs/differences_GCAT_Pandora_filtered.png")
hist(hist_subset$diff_Pandora_GCAT, main = "GCAT vs Pandora",
     xlab = "Frequency differences", ylab = "")
dev.off()
# As we can see there are some variants with very high differences

# We want to check if there are problems in the raw data but we do not find it
# check this snp with high differences
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.psam ",
              " --snp rs7627543:12631944:T:C --recode vcf",
              " --out /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/prueba"))

prueba = fread("/imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/prueba.vcf")

dim(prueba)

table(prueba[1,6000:6020]) # all 1/1

# Very big difference between Pandora and GCAT example #

freq_g[freq_g$ID %in% "2:29446721_A_G",]
#   CHROM             ID        REF   ALT     ALT_FREQS    OBS_CT
#     1:      2 2:29446721_A_G   G      A     0.024392      9976
freq_p[freq_p$ID %in% "2:29446721_A_G",]
#   CHROM             ID         REF     ALT    ALT_FREQS_PANDORA     OBS_CT_PANDORA
#     1:      2 2:29446721_A_G    G       A               0.99494               3162


# We need to create the chromosome position

#variants_library = read.xlsx("../02_Pandora_data_preparation/first_data/PandoraVariantsLibrary_20190129.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
#colnames(variants_library)[11] = "gstart"

dim(variants_library) # 47569    35

# We select the variables we are interested in
pandora_library = as.data.frame(variants_library) %>% select(Gene, 
                                                      gstart,
                                                      Ref,
                                                      Alt,
                                                      sift_pred,
                                                      polyphen2_hdiv_pred,
                                                      polyphen2_hvar_pred,
                                                      mutationassessor_pred,
                                                      mutationtaster_pred,
                                                      provean_pred)

# Delete the gstart empty variants
sum(is.na(pandora_library$gstart))
pandora_library = pandora_library[-which(is.na(pandora_library$gstart)),]

#to take the chromosome 
pandora_library$chr = as.factor(gsub("chr", "", (tstrsplit(pandora_library$gstart[grep(":",pandora_library$gstart)],split=":")[[1]])))

#to take the position
pandora_library$pos = as.integer(tstrsplit(pandora_library$gstart[grep(":",pandora_library$gstart)],split=":")[[2]])

# To create the chr_pos variable
pandora_library$chr_pos = paste0(pandora_library$chr, "_", pandora_library$pos)

# fwrite(pandora_library, "pandora_classification.txt")

## Clasification of the variants with higher difference than 0.8 between Pandora and 1000Genomes
# Now we have the data, we want to check the 
      
head(variants_library)
head(freqs)

# We want to see the classification of the higher difference variants (>0.7)
higher_differences = freqs %>% filter(diff_Pandora_GCAT > 0.7)
dim(higher_differences) # 76 VARIANTS

# Now we add the information of the classification
higher_differences = left_join(higher_differences, (variants_library %>% select(chr_pos, 
                                                                                sift_pred, 
                                                                                polyphen2_hdiv_pred,
                                                                                polyphen2_hvar_pred,
                                                                                mutationassessor_pred,
                                                                                mutationtaster_pred,
                                                                                provean_pred)))
head(higher_differences)

sum(is.na(higher_differences[,17:22])) # Almost all have not classification

# check which ones have classification
which(!is.na(higher_differences$provean_pred)) # only one
higher_differences[which(!is.na(higher_differences$provean_pred)),]

# and only have the result of mutationtaster_pred = P, which means polymorphism_automatic

# Now we want to see the same but with a filter >0.2
differences = freqs %>% filter(diff_Pandora_GCAT > 0.2)
dim(differences) # 451 VARIANTS

# Now we add the information of the classification
differences = left_join(differences, (variants_library %>% select(chr_pos, 
                                                                  sift_pred, 
                                                                  polyphen2_hdiv_pred,
                                                                  polyphen2_hvar_pred,
                                                                  mutationassessor_pred,
                                                                  mutationtaster_pred,
                                                                  provean_pred)))
head(differences)

sum(is.na(differences[,17:22])) # again almost all have not classification

# check which ones have classification
which(!is.na(differences$provean_pred)) # only three
differences[which(!is.na(differences$provean_pred)),] %>% select(chr_pos,
                                                                 sift_pred,
                                                                 polyphen2_hdiv_pred,
                                                                 polyphen2_hvar_pred,
                                                                 mutationassessor_pred,
                                                                 mutationtaster_pred,
                                                                 provean_pred)

# We have the same as in the previous filter plus:
# variant 3_10088343: Tolerant, Benign, Benign, Low (predicted as non-funcional), N (polymorphism), Neutral
# variant 3_10088343: Tolerant, Benign, Benign, neutral (predicted as non-funcional), N (polymorphism), Neutral

# All the other variants have not classification:
nrow(differences) -3 # 448

# We want to see the N of this variants with high differences

differences2 = freqs %>% filter(diff_Pandora_GCAT > 0.2) %>% select(ID,
                                                                    REF,
                                                                    ALT,
                                                                    ALT_FREQS_GCAT,
                                                                    OBS_CT_GCAT,
                                                                    Ref,
                                                                    Alt, 
                                                                    ALT_FREQS_PANDORA,
                                                                    OBS_CT_PANDORA,
                                                                    diff_Pandora_GCAT)

dim(differences2)

# Make the table
table(differences2$OBS_CT_PANDORA)

# Look the ones with < 100 observations

differences_smaller100 = differences2 %>% filter(OBS_CT_PANDORA > 100)
dim(differences_smaller100) # 88 variants observed in less than 100 individuals  
head(differences_smaller100)

# Rare variants with high observation rate
differences_high = differences2 %>% filter(OBS_CT_PANDORA > (90*3162/100))
dim(differences_high) 
head(differences_high)
table(differences_high$OBS_CT_PANDORA)
hist(differences_high$diff_Pandora_GCAT)

# Rare variants with high observation rate
differences_high2 = differences2 %>% filter(OBS_CT_PANDORA >= (3162))
dim(differences_high2) 
head(differences_high2)
table(differences_high2$OBS_CT_PANDORA)
hist(differences_high2$diff_Pandora_GCAT)


# So we filter the data to extract the variants with this difference < 0.2

# In the bim file we have all the SNPs we have previously selected and we have to delete the ones with difference > 0.2
head(bim) # We have chr_pos
dim(bim) # 13718 variants

# In differences we have all the ones with diff_Pandora_GCAT > 0.2
head(differences) # and we have also chr_pos 
dim(differences) # 451 variants

head(differences_smaller100) # IVAN

### Filter the freqs by the N_pandora

# 
dim(freqs) # 13861    15
first_final_subset = freqs %>% filter(N_Pandora > (90*3162/100))
dim(first_final_subset) # 13530    15
# With this filter we loose 331 variants

png("graphs/differences_PG_filtered.png")
hist(first_final_subset$diff_Pandora_GCAT,
     main = "Differences frequencies Pandora vs GCAT",
     xlab = "Differences") 
dev.off()

png("graphs/differences_PG_filtered.png")
hist((first_final_subset %>% filter(diff_Pandora_GCAT > 0.2))$diff_Pandora_GCAT,
     main = "Differences frequencies Pandora vs GCAT (>0.2)",
     xlab = "Differences")
dev.off()



#22_21344019 prove that pandora_nuestro is ok

selected_samples = fread("/imppc/labs/dnalab/share/PRS/02_Pandora_data_preparation/selected_merge_samples.txt")

dim(selected_samples)

aux = selected_samples %>% filter(chr %in% 22 & start %in% 21344019)

# Clasify gaps search for each gap if it is because they do not search the variable or not found it
clasify_gaps <- function(dataset){
  for(i in 7:(ncol(dataset)-1)){ # 1:1588 Sample by sample # The last column is the chr_pos
    # Take the id of the sample
    print(i)
    id_sample = colnames(dataset)[i]
    # Search the panel of this sample
    panel_sample = samples_information[samples_information$sampleId %in% 
                                         id_sample,]$version_panel
    # take as sample the genotypes of this sample
    sample = as.character(t(dataset[,i,with=F]))
    # Search for the empty genotypes on the sample: we do not know if they are NA or 0/0
    which_vacio = which(sample=="")
    # Search the empty variants 
    homos = panel[which(panel$chr_pos %in% dataset$chr_pos[which_vacio]),]
    # Search the ones with true in the correspondent panel which means this variant have been analyzed (it is 0/0)
    which_homo = as.logical(t(homos[,which(colnames(homos) %in% panel_sample),with=F]))==TRUE
    # Take the chr_pos of this 0/0 variants 
    homo_positions = homos$chr_pos[which(which_homo==T)]
    # Change this chr_pos to the 0/0
    dataset[dataset$chr_pos %in% homo_positions,i] = "0/0" 
  }
  return(dataset)
}

panel = fread("/imppc/labs/dnalab/share/PRS/04_Merge_datasets/output_data/panels_subset.txt")
aux2 = clasify_gaps(selected_samples)

aux3 = aux2 %>% filter(chr %in% 22 & start %in% 21344019)
table(as.character(aux3[,7:1587]))


interested_snps = fread("/imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/intersection_snps.txt")

head(interested_snps)

interested_snps[interested_snps$chr_pos %in% "22_21344019",]


head(first_final_subset)

# Calculate the CI 

# GCAT
first_final_subset$lowerGCAT = first_final_subset$ALT_FREQS_GCAT - 1.96 * sqrt((first_final_subset$ALT_FREQS_GCAT*(1-first_final_subset$ALT_FREQS_GCAT))/first_final_subset$N_GCAT)
first_final_subset$uperGCAT = first_final_subset$ALT_FREQS_GCAT + 1.96 * sqrt((first_final_subset$ALT_FREQS_GCAT*(1-first_final_subset$ALT_FREQS_GCAT))/first_final_subset$N_GCAT)

# Pandora
first_final_subset$lowerPandora = first_final_subset$ALT_FREQS_PANDORA - 1.96 * sqrt((first_final_subset$ALT_FREQS_PANDORA*(1-first_final_subset$ALT_FREQS_PANDORA))/first_final_subset$N_Pandora)
first_final_subset$uperPandora = first_final_subset$ALT_FREQS_PANDORA + 1.96 * sqrt((first_final_subset$ALT_FREQS_PANDORA*(1-first_final_subset$ALT_FREQS_PANDORA))/first_final_subset$N_Pandora)

?findInterval

x <- c(4,4.5)
v <- c(5, 10) # create two bins [5,10) and [10,15)
findInterval(x, v)

first_final_subset$findInterval = FALSE

for(i in 1:nrow(first_final_subset)){
  
  x = c(first_final_subset$lowerGCAT[i],first_final_subset$uperGCAT[i])
  v = c(first_final_subset$lowerPandora[i],first_final_subset$uperPandora[i])
  
  aux = findInterval(x, v)
  
  if(sum(aux %in% 1)>0){first_final_subset$findInterval[i] = TRUE}
  
  
}

table(first_final_subset$findInterval)

prova = first_final_subset[first_final_subset$findInterval,]

dim(prova)

hist(prova$diff_Pandora_GCAT)

prova %>% arrange(desc(diff_Pandora_GCAT)) %>% head()

hist(first_final_subset[676:nrow(first_final_subset),]$ALT_FREQS_PANDORA)


# Select the variables we want to discard which are the ones with differences between pandora and GCAT > 0.2
discard_variants_first_subset = first_final_subset %>% filter(diff_Pandora_GCAT >0.2) %>% select(ID)
fwrite(discard_variants_first_subset, "discard_variants_first_subset.txt")



system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.psam ",
              " --exclude discard_variants_first_subset.txt ",
              " --geno 0.1 --make-pgen ",
              " --out /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset"))



# Now we make a PCA with a filter to the rare variants (maf > 0.01)
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset.psam ",
              " --pca approx --maf 0.01 ",
              " --out /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset"))

# Read the Pca eigenvectors values 
pca_out = fread("outputs_plink/first_final_subset.eigenvec")

head(pca_out)
dim(pca_out)

# Prepare the colours
my_cols = c(rep(1,4988),rep(2,1581))

# Make the plot
png("graphs/pca2.png")
plot(pca_out$PC1,pca_out$PC2,pch=16,col=my_cols, 
     main = "Pandora vs GCAT PCA", 
     xlab= "PC1", ylab = "PC2")
legend("bottomleft",legend = c("GCAT", "Pandora"), col = c("black", "red"), 
       title = "Populations", pch = 16)
dev.off() 


# Select the variables we want to discard which are the ones with differences between pandora and GCAT > 0.05
discard_variants_first_subset = first_final_subset %>% filter(diff_Pandora_GCAT >0.01) %>% select(ID)
fwrite(discard_variants_first_subset, "discard_variants_first_subset_001.txt")



system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/prs_data.psam ",
              " --exclude discard_variants_first_subset_001.txt ",
              " --geno 0.1 --make-pgen ",
              " --out /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset_001"))



# Now we make a PCA with a filter to the rare variants (maf > 0.01)
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --pgen /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset_001.pgen ",
              " --pvar /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset_001.pvar ",
              " --psam /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset_001.psam ",
              " --pca approx --maf 0.01 ",
              " --out /imppc/labs/dnalab/share/PRS/05_GWAS/outputs_plink/first_final_subset_001"))

# Read the Pca eigenvectors values 
pca_out = fread("outputs_plink/first_final_subset_001.eigenvec")

head(pca_out)
dim(pca_out)

# Prepare the colours
my_cols = c(rep(1,4988),rep(2,1581))

# Make the plot
png("graphs/pca2_001.png")
plot(pca_out$PC1,pca_out$PC2,pch=16,col=my_cols, 
     main = "Pandora vs GCAT PCA", 
     xlab= "PC1", ylab = "PC2")
legend("bottomleft",legend = c("GCAT", "Pandora"), col = c("black", "red"), 
       title = "Populations", pch = 16)
dev.off()


## Juntar 70K SNPs de 1000G (1000Gselected.pvar) y GCAT-Pandora 3K SNPs

test = fread("things/1000Gselected.pvar",skip=188)

dim(test)
head(test)
length(unique(test$ID))

