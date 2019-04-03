# Variants imputed subset creation

In this part we are going to extract the variants of inteest from the four panels in order to have the data subset to impute with Pandora. 

We are going to select the variants of interest, which we select and save into a file *intersection_snps.txt" in the first part *01_first_data_visualization*. This file contain the information about the common variants between Pandora and the four databases. The file format is the following:

   - **snp_id**
   - rs_id
   - position
   - a0
   - a1
   - exp_freq_a1
   - **info**
   - certainty
   - type 
   - info_type0
   - concord_type0
   - r2_type0
   - **panel**
   - **chr_pos**
   
 To make this extraction we use the *--extract* function of PLINK2. We make this for the four panels selecting for each variable, from the available panels the one with the highest info. 
 
 Finally, we obtain the information only for the variants in common with Pandora in the two PLINK formats: .bed, .bim, .fam and .pgen, .psam, .pvar.  
   


