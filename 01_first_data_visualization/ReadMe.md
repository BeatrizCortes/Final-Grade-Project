## 01_First_data_visualization

The first data we receive is some information about the genes and SNPs we are going to study.

We want to check it to know more about it and have an idea of how to work with the data.

The .xlsx file we receive contain the following variables:

	- gene
	- cDNA annotation
	- Protein annotation
	- NM (Accession Number)
	- Transcrip region
	- Exon in transcript name
	- Validated cDNA
	- Validated protein
	- Effect 
	- Validated effect
	- G start (chromosome and position)
	- Ref (Referential allele)
	- Alt (Alternative allele)
	- %FP
	- inROI
	- Classification
	- Additional risk information
	- User
	- CI date (Internal Classification Data)
	- Intervar
	- Reasoning
	- Comments
	- dbSNP
	- ClinVar
	- In house frequency (frequency of all samples)
	- In house number (total number of samples containing this SNP)
	- Pop freq max (maximun population frequency)
	- N_1000G_eur (1000Genomes August 2015)
	- Exac nfe
	- Polyphen2 hdiv pred
	- Polyphen2 hvar pred
	- Mutation assessor pred
	- Mutation taster pred
	- Provean pred

We want to cross Pandora with out internal database and other open databases, to know how many common SNPs have this databases with Pandora. 

Table1. Total number of SNPs in each of the datasets we are going to use.


|                            | Total number of SNPs |
|----------------------------|----------------------|
| Pandora                    | 45,489               |
| SNParray GCAT (maf = 0.01) | 756,003              |
| Imputation (maf = 0.01)    | 15,380,370           |
| 1000Genomes (maf = 0.01)   | 14,131,756           |
| 1000Genomes (no filter)    | 82,732,419           |

Table2. Total number of SNPs in common between Pandora and the selected databases.

|                                     | Total number of SNPs |
|-------------------------------------|----------------------|
| Pandora vs SNParray (maf = 0.01)    | 496                  |
| Pandora vs imputed (maf = 0.01)     | 6,085                |
| Pandora vs 1000Genomes (maf = 0.01) | 5,182                |
| Pandora vs 1000Genomes (no filter)  | 15,340               |


This results show a small percentage of shared SNPs between Pandora and the other datasets, but this can have an explanation. 

One important point is that the other datasets are very general in comparison with Pandora. That means Pandora dataset includes only 135 genes but a lot of variants (and in general rare variants) for this genes. In the other databases we have a lot of genes but small number of variants per gene, this is easy to see in the SNPs positions, if this positions are very close, as in Pandora, that means we have a lot of close variants.

In other hand we are using filtered data: we have filters by minor allele frequency (maf) > 0.01, which means the rare SNPs where deleted, but in Pandora we have all the data without filters. Also the imputed filtered data have been filtered by info > 0.7. A histogram was made in order to see the info distribution of the common alleles. 

We also perform a PCA analysis in order to know if we loose the differentiation capacity of the SNPs between population, when we select the common SNPs between Pandora and 1000Genomes. The results shows a wide lost of differentiation capacity but still we can be able to differentiate between populations.

![1000Genomes PCA](/graphs/1000G_PCA.png)

![Common Pandora and 1000Genomes SNPs PCA](/graphs/1000G_PCA.png)


