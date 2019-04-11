######################################################################################
# 04 MERGE OF THE DATASETS (PANDORA AND GCAT)                                        #
######################################################################################

# Load libraries and set working directory
library(data.table)
library(dplyr)
library(UpSetR)

setwd("/imppc/labs/dnalab/share/PRS/04_Merge_datasets")

######################################################################################
# FIRST WE NEED TO PREPARE THE DATA                                                  #
######################################################################################

# Read the data: Pandora and GCAT_subset
pandora_merge = fread("/imppc/labs/dnalab/share/PRS/02_Pandora_data_preparation/selected_merge_samples.txt")
gcat_merge_bim = fread("/imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/outputs/final_merge/prs_gcat_imputed.bim")

# In order to cross both datasets we need a common variable, that will be chr_pos
pandora_merge$chr_pos = paste0(pandora_merge$chr, "_", pandora_merge$start)
gcat_merge_bim$V7 = paste0(gcat_merge_bim$V1, "_", gcat_merge_bim$V4)
head(gcat_merge_bim)

# Check the dimensions
dim(pandora_merge)
# 35893  1588
dim(gcat_merge_bim)
# 17252     7

# We do not want to loose information so we are going to check one by one:
#if we do not have information (snp not analyzed) (NA)
# or if they analyze the variant but do not find it (0/0)

# Read the samples information and the panels information 
samples_information = fread("../02_Pandora_data_preparation/samples_information.txt")
panels = fread("../02_Pandora_data_preparation/panels.txt")

# Change the colnames to easily work
colnames(pandora_merge) = gsub(".x","",colnames(pandora_merge))
colnames(pandora_merge) = gsub(".y","",colnames(pandora_merge))
colnames(panel)[3:6] = c("2.0","2.1","2.2","2.3")

pandora_merge[1:10,1:10]
dim(pandora_merge) #35893  1588

# Check the common variants between both
sum(gcat_merge_bim$V7 %in% pandora_merge$chr_pos) #14940
# Of the 17252 we only have 14940. 
# So we loose 2312, the reason is the two filters used: joinqualityscore and tissue

#This are the 2312 lost:
length(which(!gcat_merge_bim$V7 %in% pandora_merge$chr_pos))
lost_variants = gcat_merge_bim[which(!gcat_merge_bim$V7 %in% pandora_merge$chr_pos),]
lost_variants$V2 = NULL
head(lost_variants)

# fwrite(lost_variants, "output_data/lost_variants_2312.txt")

# SORT Pandora by chr_pos in the same order as gcat #
# In order to make the merge of the two datasets we need to have the same order on both

# We select from pandora merge only the variants that we have in gcat_merge
pandora_merge_final = pandora_merge[which(pandora_merge$chr_pos %in% gcat_merge_bim$V7),]

# There are not correspondence between the unique chr variants and the dimension
length(unique(pandora_merge_final$chr_pos)) #14940
dim(pandora_merge_final) # 15392
# This is because we have variants repeted

# Let's check:
# First we check the duplicated positions 
pandora_dup = which(duplicated(pandora_merge_final$chr_pos))
# Then obtain the chr_pos of the duplicated
pandora_dup = pandora_merge_final$chr_pos[pandora_dup]

# If we check this positions 
length(unique(pandora_dup)) # 381 variants with duplications
pandora_merge_final[pandora_merge_final$chr_pos %in% pandora_dup,][1:10,1:6]

# Remove multiallelic positions #
# this duplications are due to the fact that some positions have not single variants
pandora_merge_final = pandora_merge_final[!pandora_merge_final$chr_pos %in% pandora_dup,]

# We also filter gcat to have the same 14559 variants
gcat_merge_final = gcat_merge_bim[which(gcat_merge_bim$V7 %in% pandora_merge_final$chr_pos),]

# Check the correct dimensions
dim(gcat_merge_final) #14559     7
dim(pandora_merge_final) #14559  1588
sum(pandora_merge_final$chr_pos %in% gcat_merge_final$V7) #14559

## CREATE THE FUNCTIONS ##

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

# Create genotypes function, creates from the 1/0, 1/1, 0/0 the three columns of genotypes
create_genotypes_format <- function(dataset){
  # Also we need to obtain 3 columns for each snp with 0 0 0, 0 1 0 or 0 0 1.
  
  # Now we take the genotypes variable with the 6 first columns (the ones not genotypes)
  all_genotypes = dataset[,1:6]

  for(i in 7:(ncol(dataset)-1)){ #for each column 1588
    # We take the correspondent sample (column)
    print(i)
    sample = as.character(t(dataset[,i,with=F]))
    # Create the genotype variable with the three columns, all as 000 (all as unknown)
    geno_data = matrix(c(0,0,0),nrow = nrow(dataset),ncol = 3)  
    # Now we need to change the values of the know genotypes:
    geno_data[sample=="0/0",] = matrix(c(1,0,0), # For the 0/0 (dominant homozygote)
                                       nrow = sum(sample=="0/0"),ncol = 3,
                                       byrow=TRUE) 
    
    geno_data[sample=="1/0",] = matrix(c(0,1,0),# For the 1/0 heterozygote
                                       nrow = sum(sample=="1/0"),ncol = 3,
                                       byrow=TRUE) 
    geno_data[sample=="1/1",] = matrix(c(0,0,1), # For the 1/1 (recesive homozygote)
                                       nrow = sum(sample=="1/1"),ncol = 3,
                                       byrow=TRUE) 
    all_genotypes = cbind(all_genotypes,geno_data) # Now bind the new 3 column created to the 6 first variables
  }
  return(all_genotypes)
}




# FOR PANDORA_MERGE without the filter of the data
start = Sys.time()
pandora_allmerge = clasify_gaps(pandora_merge)
end = Sys.time()
time= end - start
time # Time difference of 7.082442 mins

# Also we need to obtain 3 columns for each snp with 0 0 0, 0 1 0 or 0 0 1.

# Now we take the pandora_allgenotypes variable with the 6 first columns (the ones not genotypes)

start = Sys.time()
pandora_allgenotypes = create_genotypes_format(pandora_allmerge)
end = Sys.time()
time= end - start
time # Time difference of 13.98865 mins

# Check
pandora_allgenotypes[1:10,1:10]
dim(pandora_allgenotypes) # 35893  4749:  1588 -7  = 1581 individuals. 1581*3 (3 for each sample) +6 (the variables)

# Write
# fwrite(pandora_allgenotypes, "output_data/pandora_allmerge_genotypes_35893.txt")
# pandora_allgenotypes = fread("output_data/pandora_allmerge_genotypes_35893.txt")

# Now REPEAT with the SELECTED SUBSET
start = Sys.time()
pandora_merge_subset = clasify_gaps(pandora_merge_final)
end = Sys.time()
time= end - start
time # Time difference of 2.558158 mins

dim(pandora_merge_final) #13718  1588
table(pandora_merge_final[,1587])
#0/0   1/0   1/1 
#349 13001   657   552 

# Also we need to obtain 3 columns for each snp with 0 0 0, 0 1 0 or 0 0 1.

# Now we take the pandora_genotypes variable with the 6 first columns (the ones not genotypes)

start = Sys.time()
pandora_genotypes = create_genotypes_format(pandora_merge_subset)
end = Sys.time()
time= end - start
time # Time difference of 6.701114 mins

# Check
pandora_genotypes[1:10,1:10]
dim(pandora_genotypes) # 13718  4749:  1588 -7  = 1581 individuals. 1581*3 (3 for each sample) +6 (the variables)

# write
# fwrite(pandora_genotypes, "output_data/pandora_subset_genotypes_13718.txt")
# pandora_genotypes = fread("output_data/pandora_subset_genotypes_13718.txt")

# Read genotype probabilities #
# Read the .gen file
my_gen = fread("/imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/outputs/final_merge/my_files.gen")
my_gen[1:5,1:6]
dim(my_gen) #17252 14969
# Again it have all the variants (17252) without the filters

# So, we filter the data:
my_gen = my_gen[my_gen$V2 %in% gcat_merge_final$V2,] # 2693 filtered
dim(my_gen) # 14559 14969
my_gen[1:5,1:6]

# Add the new common variable
my_gen$chr_pos = paste0(my_gen$V1, "_", my_gen$V3)

# We check if the order of the variants is the same
identical(pandora_merge_final$chr_pos,my_gen$chr_pos) # FALSE

# ORDER BY CHR_POS #
# We order all by chr_pos
my_gen = my_gen[order(my_gen$chr_pos),]
gcat_merge_final = gcat_merge_final[order(gcat_merge_final$V7),]
pandora_merge_final = pandora_merge_final[order(pandora_merge_final$chr_pos),]

# Check if they are true
identical(pandora_merge_final$chr_pos,my_gen$chr_pos) # TRUE
identical(gcat_merge_final$V7,my_gen$chr_pos) # TRUE

# CHECK REFERENCE ALLELES #
# We need to check if the reference and the alternative allele is the same
# If it is not the same we need to change the bases and the genotype.

# It is very important
# we need to make the merge with the same reference and alternative allele for all the samples
pandora_merge_final[1:10,1:10]
my_gen[1:10,1:10]

# We check if it is the same or the inverse order 
identical(pandora_merge_final$refBase, my_gen$V4) # FALSE. ref != ref
identical(pandora_merge_final$refBase, my_gen$V5) # FALSE. ref != alt

# So, we need to check all variants one by one and determine if it is equal, the inverse or it is wrong

# Identify changes
# we identify the variants that have inversions: ref1 = alt2 and alt1 = ref2
snps_change = pandora_merge_final[which((pandora_merge_final$refBase == my_gen$V5) & 
                (pandora_merge_final$altBase == my_gen$V4)),]
dim(snps_change) # We have 3 inverted variants.

# First we change the genotypes 
snps_change[snps_change=="1/1"] = "a/a"
snps_change[snps_change=="0/0"] = "1/1"
snps_change[snps_change=="a/a"] = "0/0"

snps_change[,1:10]

# Then change the ref and alt bases
a = snps_change$refBase
snps_change$refBase = snps_change$altBase
snps_change$altBase = a

# Now we need to make another subset only with the rigth ones
snps_iden = pandora_merge_final[which((pandora_merge_final$refBase == my_gen$V4) & 
                                          (pandora_merge_final$altBase == my_gen$V5)),]
dim(snps_iden) # 13715
nrow(snps_iden) + nrow(snps_change) # nos quedamos solo con REF=REF, ALT=ALT y  REF=ALT, ALT=REF

dim(pandora_merge_final) # 14559  1588
dim(snps_change) # 3 1588
dim(snps_iden) # 13715  1588
14559 - 13715 # 844

# Now we merge both the identical and oposite datasets and order it again 
pandora_merge_final = rbind(snps_iden,snps_change)
pandora_merge_final = pandora_merge_final[order(pandora_merge_final$chr_pos),]
my_gen = my_gen[my_gen$chr_pos %in% pandora_merge_final$chr_pos,]
my_gen = my_gen[order(my_gen$chr_pos),]

# Check if they are identical: the same variant order
identical(my_gen$V4,pandora_merge_final$refBase) # TRUE
identical(my_gen$V5,pandora_merge_final$altBase) # TRUE
identical(my_gen$chr_pos,pandora_merge_final$chr_pos) # TRUE

# Prepare pandora_merge_final
# In pandora_merge_final we have 1/1, 1/0, 0/0 in one column
# We need three columns with  0 0 1, 0 1 0 or 1 0 0.

# Take care! We same some missing values as "" 
pandora_merge_final[1:5,1:10]

# Create the new data frame: pandora_final
pandora_final = pandora_merge_final[,1:6]
pandora_genotypes = pandora_merge_final[,7:ncol(pandora_merge_final)]
pandora_genotypes = pandora_genotypes %>% select(-chr_pos)


# Panels:
panels$panel20 = ifelse(panels$panel20==T,1,0)
panels$panel21 = ifelse(panels$panel21==T,1,0)
panels$panel22 = ifelse(panels$panel22==T,1,0)
panels$panel23 = ifelse(panels$panel23==T,1,0)


png("graphs/all_snps_by_panel.png", width = 16,height = 8,res = 300,units = "in")
upset(panels,nsets = 4,sets = c("panel20","panel21","panel22","panel23"),
      order.by = c("freq"), decreasing = c(TRUE))
dev.off()

# Now we want to make a subset of the SNPs, only with the 13718 SNPs we are interested in 

dim(pandora_final) # 13718     6
head(pandora_final)
# We take in the variable snps_subset all the variants on pandora_final
snps_subset = paste0(pandora_final$chr, "_", pandora_final$start)
head(snps_subset)
length(snps_subset)#13718

# Now we filter the data: panels_subset
panels_subset = panels %>% filter(panels$chr_pos %in% snps_subset)
dim(panels_subset) # 13747     6

png("graphs/subset_snps_by_panel.png", width = 16,height = 8,res = 300,units = "in")
upset(panels_subset,nsets = 4,sets = c("panel20","panel21","panel22","panel23"),
      order.by = c("freq"), decreasing = c(TRUE))
dev.off()

#fwrite(panels_subset, "output_data/panels_subset.txt")

##11753 SNPs are common between the four panels:
13718- 11753 # 1965

# we want to know if we loose information of the genes if we delete this 1965 SNPs
subset_4panels = data.frame("Gene" = "", "chr_pos" = "") # create the new data frame
for (r in 1:nrow(panels_subset)){ # For each row
  snps_actual = panels_subset[r,] # take the correspondent snp
  if (rowSums(snps_actual[,3:6]) == 4){ # if the snps is in the four panels 
    subset_4panels = rbind(subset_4panels, snps_actual[,1:2]) # add the row to the data.frame
  }
}

length(unique(subset_4panels$Gene)) # 116 genes in common (in panels_subset we have 137)
table(subset_4panels$Gene) # we loose 21 genes, 

######################################################################################
# NOW WE CAN MAKE THE MERGE                                                          #
######################################################################################

# We want to add on the rigth of the my_gen the genotypes of Pandora
# We need to take into account:
#     dim my_gen = 13718 x 14970 (14970-5variables - 1 (chr_pos)) = 14964 SAMPLES
#     dim pandora_merge_final  = 13718 x 1588 (1588 - 6 variables - 1 (chr_pos)) = 1581 SAMPLES
#     in both we have the same number of SNPs (13718) but 14964 samples vs 1581

# check the two datasets have the same SNPs order
identical(my_gen$chr_pos,pandora_merge_final$chr_pos) # TRUE

# Prepare the two datasets:
# delete the chr_pos from my_gen
my_gen$chr_pos = NULL

# Delete all the variables of pandora 
pandora_merge_final = pandora_merge_final[,7:ncol(pandora_merge_final), with=F]
pandora_merge_final$chr_pos = NULL

# Now we can make the merge:
my_final_genotypes = cbind(my_gen, pandora_merge_final)

#Check
dim(my_final_genotypes) # 13718 x 16550

my_final_genotypes[1:10,1:10]

# fwrite(my_final_genotypes, "output_data/my_final_genotypes.gen")
# fwrite(my_final_genotypes, "output_data/my_final_genotypes.txt")

#To make the pgen file
system(paste0(" /imppc/labs/dnalab/share/PRS/03_GCAT_subset_variants_imputed/plink2 ",
              " --gen /imppc/labs/dnalab/share/PRS/04_Merge_datasets/output_data/my_final_genotypes.gen ",
              " --sample /imppc/labs/dnalab/share/data_storage/GCATCore_imputed_raw_august_18/GCAT_imputation_August_2018/common/mixed/Chr_1/mixed_phasing_chr_1.sample ",
              " --make-pgen ",
              " --sort-vars ",
              " --out /imppc/labs/dnalab/share/PRS/04_Merge_datasets/prs_data"))
# Error: Line 1 of .gen file has fewer tokens than expected.

#To make the bed file
system(paste0("plink2 ",
              " --gen /imppc/labs/dnalab/share/PRS/04_Merge_datasets/output_data/my_final_genotypes.gen ",
              " --sample /imppc/labs/dnalab/share/data_storage/GCATCore_imputed_raw_august_18/GCAT_imputation_August_2018/common/mixed/Chr_1/mixed_phasing_chr_1.sample ",
              " --make-bed ",
              " --out /imppc/labs/dnalab/share/PRS/04_Merge_datasets/prs_data"))
# Error: Invalid chromosome code
table(my_final_genotypes$V1)

