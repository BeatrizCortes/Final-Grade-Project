################################################################################
#      PANDORA SECOND STEP: DATA MANAGEMENT AND VISUALIZATION                  #
################################################################################

#### Libraries ####
library(dplyr)
library(data.table)

#### Working directory ####
setwd("/imppc/labs/dnalab/share/PRS/02_Pandora_data_preparation")

##### Load the data ####

### all.samples contain the information about each sample 
# it contain control and cases
load("./AllSamples.RData")
head(all.samples)

# the important variants are:
#     - SampleName/individualProgenyId
#     - sampleId
#     - version_panel (version del panel)
#     - tissue (we are only interested in blood)
names(all.samples)[16] <- "version_panel"

### control.samples contain the control samples
load("control.samples.RData")
head(control.samples)

# the sampleId of the control samples
# It is important to label the samples in order to know which are case/control

### all.vars contain the variants index by analysisId
# It is a list where each element of the list is one sample/patient
load("all.vars.RData")

# genes.in.design.versions to know the genes included in each design
load("genes.in.design.versions.RData")
# it is a list of 5 elements one for each design

##############################################################################
# DATA MANAGEMENT                                                            #
##############################################################################

# We need a loop with the same number of iterations as the number of samples = 2530
# In each iteration we will select only the ones with "sangre" as tissue
table(all.samples$tissue) # sangre: 2142

#We check if they have the same order
identical(names(all.vars), all.samples$analysisId) # no same order 

# Create the vars variable with all the vars names 
vars = names(all.vars) 
# With this we take the analysisId of each list inside all.vars

# Auxiliar function: We need to change the freq alelica to the genotypes
obtain_genotypes <- function(sample1){
  sample1$genotype = sample(0, nrow(sample1), replace = T)
  for (j in 1:nrow(sample1)){
    frequency = sample1[j,]$freq
    if (frequency >= 80){ #if the frequency is > 80 sustitute by genotype 1/1 ALt homozygote
      sample1[j,]$genotype = "1/1"
    }
    if (frequency <= 20){#if the frequenct is < 20 sustitute by genotype 0/0 REF homozygote
      sample1[j,]$genotype = "0/0"
    }
    if (frequency < 80 && frequency > 20){#else it is heterocygote, so genotype 1/0
      sample1[j,]$genotype = "1/0"
    }
  }
  return(sample1)
}


# Principal loop: merge the data
#Divide the vars, because if not the loop takes a lot of time
vars1 = vars[1:500]
vars2 = vars[501:1000]
vars3 = vars[1001:1500]
vars4 = vars[1501:2000]
vars5 = vars[2001:2500]
vars6 = vars[2501:2628]


# 1st STEP # (using vars1)
merge_samples = NULL
c = 0
start = Sys.time()
for(i in vars1){
  # Tissue condition from all.samples
  # Check if the row in all.samples with the same analysisId have as tissue Sangre
  # Also check the analysisId is on all.samples
  if ((all.samples[((all.samples$analysisId) == i),]$tissue == "Sangre") && (!is.na((all.samples[((all.samples$analysisId) == i),]$tissue == "Sangre")))){
    # All one to the counter and print it no know in which iteration we are
    c = c + 1
    print(c)
    # Now i have the samle from all.vars, the variable actual_sample
    actual_sample = all.vars[[c]]
    if (length(all.vars[[c]]) > 0){
      # But we want to filter for jointQualityScore < 4
      actual_sample = actual_sample %>% filter(jointQualityScore<4)
      # Obtain the genotypes (allele frequence) 
      actual_sample = obtain_genotypes(actual_sample)
      #With this function we obtain the genotypes
      colnames(actual_sample)[14] = unique(actual_sample$sampleId)
      if (c == 1){ #In the case of the first one we create the sample variable
        merge_samples = actual_sample%>% select(-id,-sampleId,-analysisId,-genomeVersionId,-jointQualityScore, -freq, -depth)
      }
      else{ # make the merge
        merge_samples = full_join(merge_samples,
                                  actual_sample %>% select(-id,-sampleId,-analysisId,-genomeVersionId,-jointQualityScore, -freq, -depth))
      }
    }
  }
}
end = Sys.time()
total_time1 = end - start # we calculate the total time it takes to perform the 500 iterations
# Time difference of 7.143542 mins
# dim(merge_samples) 13426   451

# Save the file
fwrite(merge_samples,"things/pandora_data_merge1.txt")

# Read the data 
merge_samples1 = fread("things/pandora_data_merge1.txt")
dim(merge_samples1) # 3426   451
head(merge_samples1)


merge_samples = NULL
##############################
 
c #446
threshold = c+1

# 2nd STEP #
start = Sys.time()
for(i in vars2){
  # tissue condition from all samples
  if ((all.samples[((all.samples$analysisId) == i),]$tissue == "Sangre") && (!is.na((all.samples[((all.samples$analysisId) == i),]$tissue == "Sangre")))){
    c = c + 1
    print(c)
    # Now i have the samle from all.vars, the variable actual_sample
    actual_sample = all.vars[[c]]
    if (length(all.vars[[c]]) > 0){
      # But we want to filter for jointQualityScore < 4
      actual_sample = actual_sample %>% filter(jointQualityScore<4)
      # Obtain the genotypes (allele frequence) 
      actual_sample = obtain_genotypes(actual_sample)
      #With this function we obtain the genotypes
      colnames(actual_sample)[14] = unique(actual_sample$sampleId)
      if (c == threshold){ #In the case of the first one we create the sample variable
        merge_samples = actual_sample%>% select(-id,-sampleId,-analysisId,-genomeVersionId,-jointQualityScore, -freq, -depth)
      }
      else{ # make the merge
        merge_samples = full_join(merge_samples,
                                  actual_sample %>% select(-id,-sampleId,-analysisId,-genomeVersionId,-jointQualityScore, -freq, -depth))
      }
    }
  }
}
end = Sys.time()
total_time2 = end - start
#Time difference of 7.181696 mins

# Save the file
fwrite(merge_samples,"things/pandora_data_merge2.txt")

# Read the data 
merge_samples2 = fread("things/pandora_data_merge2.txt")
dim(merge_samples2)

merge_samples = NULL
##############################

c #897
threshold = c+1

# 3rd STEP #
start = Sys.time()
for(i in vars3){
  # tissue condition from all samples
  if ((all.samples[((all.samples$analysisId) == i),]$tissue == "Sangre") && (!is.na((all.samples[((all.samples$analysisId) == i),]$tissue == "Sangre")))){
    c = c + 1
    print(c)
    # Now i have the samle from all.vars, the variable actual_sample
    actual_sample = all.vars[[c]]
    if (length(all.vars[[c]]) > 0){
      # But we want to filter for jointQualityScore < 4
      actual_sample = actual_sample %>% filter(jointQualityScore<4)
      # Obtain the genotypes (allele frequence) 
      actual_sample = obtain_genotypes(actual_sample)
      #With this function we obtain the genotypes
      colnames(actual_sample)[14] = unique(actual_sample$sampleId)
      if (c == threshold){ #In the case of the first one we create the sample variable
        merge_samples = actual_sample%>% select(-id,-sampleId,-analysisId,-genomeVersionId,-jointQualityScore, -freq, -depth)
      }
      else{ # make the merge
        merge_samples = full_join(merge_samples,
                                  actual_sample %>% select(-id,-sampleId,-analysisId,-genomeVersionId,-jointQualityScore, -freq, -depth))
      }
    }
  }
}
end = Sys.time()
total_time3 = end - start
# Time difference of 4.954592 mins

# Save the file
fwrite(merge_samples,"things/pandora_data_merge3.txt")

# Read the data 
merge_samples3 = fread("things/pandora_data_merge3.txt")
dim(merge_samples3) #13040   321

merge_samples = NULL
##############################

c # 1212
threshold = c+1

# 4th STEP #
start = Sys.time()
for(i in vars4){
  # tissue condition from all samples
  if ((all.samples[((all.samples$analysisId) == i),]$tissue == "Sangre") && (!is.na((all.samples[((all.samples$analysisId) == i),]$tissue == "Sangre")))){
    c = c + 1
    print(c)
    # Now i have the samle from all.vars, the variable actual_sample
    actual_sample = all.vars[[c]]
    if (length(all.vars[[c]]) > 0){
      # But we want to filter for jointQualityScore < 4
      actual_sample = actual_sample %>% filter(jointQualityScore<4)
      # Obtain the genotypes (allele frequence) 
      actual_sample = obtain_genotypes(actual_sample)
      #With this function we obtain the genotypes
      colnames(actual_sample)[14] = unique(actual_sample$sampleId)
      if (c == threshold){ #In the case of the first one we create the sample variable
        merge_samples = actual_sample%>% select(-id,-sampleId,-analysisId,-genomeVersionId,-jointQualityScore, -freq, -depth)
      }
      else{ # make the merge
        merge_samples = full_join(merge_samples,
                                  actual_sample %>% select(-id,-sampleId,-analysisId,-genomeVersionId,-jointQualityScore, -freq, -depth))
      }
    }
  }
}
end = Sys.time()
total_time4 = end - start
# Time difference of 6.912758 mins

# Save the file
fwrite(merge_samples,"things/pandora_data_merge4.txt")

# Read the data 
merge_samples4 = fread("things/pandora_data_merge4.txt")
dim(merge_samples4)

merge_samples = NULL
##############################

# 5th STEP #

c # 1671
threshold = c+1

start = Sys.time()
for(i in vars5){
  # tissue condition from all samples
  if ((all.samples[((all.samples$analysisId) == i),]$tissue == "Sangre") && (!is.na((all.samples[((all.samples$analysisId) == i),]$tissue == "Sangre")))){
    c = c + 1
    print(c)
    # Now i have the samle from all.vars, the variable actual_sample
    actual_sample = all.vars[[c]]
    if (length(all.vars[[c]]) > 0){
      # But we want to filter for jointQualityScore < 4
      actual_sample = actual_sample %>% filter(jointQualityScore<4)
      # Obtain the genotypes (allele frequence) 
      actual_sample = obtain_genotypes(actual_sample)
      #With this function we obtain the genotypes
      colnames(actual_sample)[14] = unique(actual_sample$sampleId)
      if (c == threshold){ #In the case of the first one we create the sample variable
        merge_samples = actual_sample%>% select(-id,-sampleId,-analysisId,-genomeVersionId,-jointQualityScore, -freq, -depth)
      }
      else{ # make the merge
        merge_samples = full_join(merge_samples,
                                  actual_sample %>% select(-id,-sampleId,-analysisId,-genomeVersionId,-jointQualityScore, -freq, -depth))
      }
    }
  }
}
end = Sys.time()
total_time5 = end - start
# Time difference of 30.21086 mins

# Save the file
fwrite(merge_samples,"things/pandora_data_merge5.txt")

# Read the data 
merge_samples5 = fread("things/pandora_data_merge5.txt")
dim(merge_samples5)

merge_samples = NULL
##############################

c # 
threshold = c+1

# 6th STEP #
start = Sys.time()
for(i in vars6){
  # tissue condition from all samples
  if ((all.samples[((all.samples$analysisId) == i),]$tissue == "Sangre") && (!is.na((all.samples[((all.samples$analysisId) == i),]$tissue == "Sangre")))){
    c = c + 1
    print(c)
    # Now i have the samle from all.vars, the variable actual_sample
    actual_sample = all.vars[[c]]
    if (length(all.vars[[c]]) > 0){
      # But we want to filter for jointQualityScore < 4
      actual_sample = actual_sample %>% filter(jointQualityScore<4)
      # Obtain the genotypes (allele frequence) 
      actual_sample = obtain_genotypes(actual_sample)
      #With this function we obtain the genotypes
      colnames(actual_sample)[14] = unique(actual_sample$sampleId)
      if (c == threshold){ #In the case of the first one we create the sample variable
        merge_samples = actual_sample%>% select(-id,-sampleId,-analysisId,-genomeVersionId,-jointQualityScore, -freq, -depth)
      }
      else{ # make the merge
        merge_samples = full_join(merge_samples,
                                  actual_sample %>% select(-id,-sampleId,-analysisId,-genomeVersionId,-jointQualityScore, -freq, -depth))
      }
    }
  }
}
end = Sys.time()
total_time6 = end - start
#Time difference of 1.549429 mins


# Save the file
fwrite(merge_samples,"things/pandora_data_merge6.txt")

# Read the data 
merge_samples6 = fread("things/pandora_data_merge6.txt")
dim(merge_samples6)

merge_samples = NULL

##############################

###########FINAL MERGE##############

final_merge_samples = full_join(merge_samples1,merge_samples2)

final_merge_samples = full_join(final_merge_samples,merge_samples3)

final_merge_samples = full_join(final_merge_samples,merge_samples4)

final_merge_samples = full_join(final_merge_samples,merge_samples5)

final_merge_samples = full_join(final_merge_samples,merge_samples6)


# Save the file
fwrite(final_merge_samples,"pandora_data_finalmerge.txt")

# Read the data 
final_merge_samples = fread("pandora_data_finalmerge.txt")

## PROBLEM empty elements of the list

length(all.vars)

empty = c()
j = 230
for (j in 1:length(all.vars)){
  print(j)
  if (length(all.vars[[j]]) < 1){
    empty = c(empty, j)
  }
}

all.vars[[230]]
all.vars[[1518]]

#total merge = 2246 (2240 ind + 6 variables) 
# Sangre = 2242 
# all.vars = 2628
# all.samples = 2628

# The two missing individuals are the two empty list individuals of all.vars

# We can see there are a lot of missing values in some samples
table(final_merge_samples$`2542`) #34728 missings

# The reason of that can be the fact that we have 4 different panels:
# for some samples we have not analyced the SNPs correspondent to some genes

# Let's check for which panel this sample is

head(all.samples)

all.samples[which(all.samples$sampleId == 2542),]$version_panel # 2.1

# this version have only 126 genes genes, when the last version have 136

# Let's check if a sample from the last version have missings

head(all.samples[which(all.samples$version_panel == "2.3"),]$sampleId) 
# [1] 6664 6679 6683 6697 6696 6699

table(final_merge_samples$`6664`) # Also have missings
#         0/0   1/0   1/1 
#34543     1   818   531 

which(all.samples$sampleId == "6664") #1919

length(all.vars[[1919]]$freq) # 1430 we have this frequences 
length(which(final_merge_samples$`6664` == ""))
#[1] 34543

# If we sum the other

length(which(final_merge_samples$`6664` != ""))
# [1] 1350

# Lets check: repeat the process step by step with this sample that we know why we loose genotypes
c = which(all.samples$sampleId == "6664")
actual_sample = all.vars[[c]]
dim(actual_sample) #1430
actual_sample = actual_sample %>% filter(jointQualityScore<4)
dim(actual_sample) #1350 

# We loose this genotypes fue to the jointQualityScore, so all is OK. 
# We have a lot of missings but this is not a problem of the data management
# we use all the freqs they give us in Pandora.

# Missings are due to the Pandora database

### all.genes ### 

# the file is a list of 5 elements:
#     1st element correspondent to panel_version 2.0 has 122 genes
#     2nd element correspondent to panel_version 2.1 has 126 genes
#     3rd element correspondent to panel_version 2.2 has 135 genes
#     4th element correspondent to panel_version 2.3 has 136 genes
#     5th element contain 3 genes: APOB, LDLR and PCSK9

head(all.genes[[1]]) # The second column is the symbol

genes_panel20 = all.genes[[1]]$symbol
genes_panel21 = all.genes[[2]]$symbol
genes_panel22 = all.genes[[3]]$symbol
genes_panel23 = all.genes[[4]]$symbol

# 1st panel improvement 
# With this we can see the genes added in the first panel improvement
first_improvement = setdiff(genes_panel21, genes_panel20)
# "LZTR1", "GRB2",  "RASA2", "PIK3CA"

# 2nd panel improvement 
second_improvement = setdiff(genes_panel22, genes_panel21)
# "A2ML1" ,"FAN1", "NTHL1", "POT1", "RAD51B", "RIT1", "RRAS", "SMARCA4", "SOS2"

# 3th improvement 
third_improvement = setdiff(genes_panel23, genes_panel22)
# "EGFR", "ERBB2", "CDKN2B", "PPP1CB","SDHA" 

# Filter the samples from the selected tissue and select the interesting variables
samples_information = all.samples %>% filter(tissue == "Sangre") %>% select(analysisId, sampleId, version_panel, GeneCascadesName)

head(samples_information)

# Now we want to add a new column to differenciate between control and case
# First we create the column and call all Case
samples_information$Type = sample("Case", nrow(samples_information), replace = T)

# Now change the controls
control = control.samples$sampleId
length(control) #239

for (c in control){
  samples_information$Type[samples_information$sampleId == c] <- "Control"
}

table(samples_information$Type) # from the 239 controls only 32 have "Sangre" as tissue

## Now we have to delete the missing information individuals
all.vars[230] # 3031 this is the analysisId
all.vars[1518] # 4594

# Check if we have in out dataset this missing samples
dim(samples_information) # 2242    5
which(samples_information$analysisId == 3031) #we need to delete the ind 205
which(samples_information$analysisId == 4594) #we alrealdy do not have this individual (we filtered it before)

# Delete the missing samples
samples_information = samples_information[-which(samples_information$analysisId == 3031),]

dim(samples_information) # We loose one individual

# Save the file
fwrite(samples_information,"samples_information.txt")

# Read the data 
samples_information = fread("samples_information.txt")


## Pandora Genes information ##

# These are the genes contained in Pandora
genes <- fread("../01_first_data_visualization/data/pandora_genes.txt")

# We only have the names and we want:
# the start_position, the end_position, the chromosome name and the strand
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)

# Create the empty data frame
genes_info = data.frame(hgnc_symbol = character(),
                        ensembl_gene_id = character(),
                        chromosome_name = integer(),
                        start_position = integer(),
                        end_position = integer(), 
                        strand = integer())

# Make a loop to search for the information of each gene and add it to the dataframe
for (gene in genes){
  genes_info = rbind(genes_info, getBM(attributes=c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"), # we indicate the information we want to retrieve
        filters=c("hgnc_symbol"), # We want to filter only the genes contained on Pandora
        values=list(gene),
        mart=ensembl))
}
length(unique(genes_info$hgnc_symbol))

# We delete the repeated genes
genes_info = genes_info[-grep("CHR", genes_info$chromosome_name),]

# Now change the strand
head(genes_info)
genes_info$strand = gsub(-1, "-", genes_info$strand)
genes_info$strand = gsub(1, "+", genes_info$strand)
head(genes_info)

# Save the file
fwrite(genes_info,"genes_information.txt")

genes_info = fread("genes_information.txt")

# We have finally 139 unique genes
length(unique(genes_info$ensembl_gene_id)) #139

## PROBLEM ## 
length(genes_panel23) # 136

# If we have 136 genes in the biggest panel 
# why we have 139 unique genes in pandora_genes?

# We know these are the new genes, but we want to know if they delete any gene.

# Prepare the data for the left join
genes_panel20 = as.data.frame(genes_panel20)
colnames(genes_panel20) = "genes"
genes_panel20$panel = sample("2.0", length(genes_panel20), replace = TRUE)

genes_panel21 = as.data.frame(genes_panel21)
colnames(genes_panel21) = "genes"
genes_panel21$panel = sample("2.1", length(genes_panel21), replace = TRUE)

genes_panel22 = as.data.frame(genes_panel22)
colnames(genes_panel22) = "genes"
genes_panel22$panel = sample("2.2", length(genes_panel22), replace = TRUE)

genes_panel23 = as.data.frame(genes_panel23)
colnames(genes_panel23) = "genes"
genes_panel23$panel = sample("2.3", length(genes_panel23), replace = TRUE)


all_genes_by_panel = NULL
# We have genes in panel 2.2 that we do not have in 2.3
not_in_2.3 = setdiff(genes_panel22$genes, genes_panel23$genes)
# We need to add them to genes_panel23
new = data.frame("genes" = not_in_2.3, "panel" = NA)
genes_panel23 = rbind(genes_panel23, new)

all_genes_by_panel = left_join(genes_panel23, genes_panel22, by = "genes")
all_genes_by_panel = left_join(all_genes_by_panel, genes_panel21, by = "genes")
genes_panel20$genes = as.character(genes_panel20$genes)
all_genes_by_panel = left_join(all_genes_by_panel, genes_panel20, by = "genes")

colnames(all_genes_by_panel) = c("genes", "panel3", "panel2", "panel1", "panel0")

dim(all_genes_by_panel)
dim(genes_panel20) # 122
dim(genes_panel21) # 126 "LZTR1", "GRB2",  "RASA2", "PIK3CA"
dim(genes_panel22) # 135 "A2ML1" ,"FAN1", "NTHL1", "POT1", "RAD51B", "RIT1", "RRAS", "SMARCA4", "SOS2"
dim(genes_panel23) # 136+4 "EGFR", "ERBB2", "CDKN2B", "PPP1CB","SDHA" 

length(all_genes_by_panel$panel0) - sum(is.na(all_genes_by_panel$panel0))
length(all_genes_by_panel$panel1) - sum(is.na(all_genes_by_panel$panel1))
length(all_genes_by_panel$panel2) - sum(is.na(all_genes_by_panel$panel2))
length(all_genes_by_panel$panel3) - sum(is.na(all_genes_by_panel$panel3))

all_genes_by_panel

fwrite(all_genes_by_panel, "panels.txt")

all_genes_by_panel = fread("panels.txt")

##############################################################################
# VISUALIZATION                                                              #
##############################################################################

## Histogram chromosomes ## 
# We want to know how many genes there are for each chromosome:
hist_genes = gsub("X", 23, genes_info$chromosome_name)
hist(as.integer(hist_genes))

png("graphs/genes_per_chr.png")
hist(as.integer(hist_genes), main = "Genes per chromosome number", 
     xlab = "chromosome", ylab = "Numer of genes", breaks = 23)
abline(h = 5, col = "cyan3", lwd=5, lty=3)
abline(h = 10, col = "cyan3", lwd=5, lty=3)
abline(h = 15, col = "cyan3", lwd=5, lty=3)
abline(h = 20, col = "cyan3", lwd=5, lty=3)
dev.off()

## KaryoplotR ##

# Install
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#sBiocManager::install("karyoploteR")

# Citation
#Bernat Gel & Eduard Serra. (2017). karyoploteR: an R/Bioconductor package to plot customizable genomes displaying arbitrary data. Bioinformatics, 31–33. doi:10.1093/bioinformatics/btx346

library("karyoploteR", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")

# We need to divide the information in the different chromosomes 
all_chr = levels(as.factor(genes_info$chromosome_name))

# Loop to make the plots

# First we need to create the vector with the chromosomes names:
chromosomes = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                "22", "X")

for (c in chromosomes){
  chr_actual = genes_info[genes_info$chromosome_name == c,]
  png(paste0("graphs/Chr", c ,"_genes.png"), width = 1600,height = 800,res = 300)
  genes = data.frame(chr=rep(paste0("chr", c), nrow(chr_actual)), pos = chr_actual$start_position, labels= chr_actual$hgnc_symbol)
  kp <- plotKaryotype(chr = paste0("chr", c))
  kpAddBaseNumbers(kp)
  kpPlotMarkers(kp, chr=genes$chr, x=genes$pos, labels=genes$labels)
  dev.off()
}