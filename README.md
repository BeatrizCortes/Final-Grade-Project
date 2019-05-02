# Final Grade Project
A repository where I will include all the scripts used during the development of the final grade project. We divide the work done in the following steps:

<br/>

1. **First data visualization**
   *In this section we check the primary information we had about the variants contained in Pandora: maf, common variants with the available datasets, info of the common variants, differenciation of populations capacity of this variants...* 
   
   <br/>
   
2. **Pandora data preparation**
   *We extract the important information from the four tables given to us and preparation of the information of samples, SNPs, genotypes...*
   
   *We check the number of genes by chromosome and the position of the genes in the different chromosomes.*
      
   <br/>
   
3. **GCAT subset variants imputed**
   *We extract the variants of interest determined on the first data visualization from the four panels, taking the ones with the higher info.*
   
   <br/>
   
4. **Merge datasets**
   *We prepare the data from Pandora and from the Subset previously obtained in order to merge them into a unique file with the information of all the selected variants (13718) and all the samples (16545).*

   <br/>

5. **GWAS Preparation** 
  *In this section we will check the data to see if it is correct. We check GCAT and Pandora data and check also the differences between the frequences of the selected SNPs.
  
  6. **GWAS**
  *We divide the data in two sets: training 70% of the samples and the test 30% of the samples. We perform the GWS for the four groups and the training set and select the correlated SNPs and check the number of genes it contain.
  
  
  7. **PRS**
  *We want to calculate predictive models with the selected SNPs as correlated for each of the cancer groups. We obtain boxplots, densityplots, contingency tables, AUC, sensitivity and specificity for each group.*
