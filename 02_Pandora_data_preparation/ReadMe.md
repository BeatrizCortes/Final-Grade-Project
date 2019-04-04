# Pandora data preparation

## Files information

We obtain the data divided into 4 different files: all.samples, all.vars, control.samples and all.genes.

We want to check the files in order to understand the information it contain and then we want to search and prepare the data we are interested in.

#### all.samples
This file contain the information about each sample with the following variables, the most important are the bolds:
  - name                  
  - laboratoryId        
  - **individualProgenyId** 
  - **individualId**         
  - SampleName           
  - SampleType           
  - **tissue**
  - **sampleId**            
  - sampleCaptureId
  - sampleDataId         
  - AnalysisName
  - PipelineName         
  - version
  - pipelineId           
  - DesignVersionLabel
  - **version_panel**              
  - DesignVersionName
  - designVersionId      
  - runId
  - RunName
  - RunDate
  - readLength
  - RunStatus
  - instrumentId
  - InstrumentName
  - InstrumentType
  - InstrumentInstitution
  - AnalysisStatus
  - AnalysisDate
  - analysisId
  - **GeneCascadesName**
  - geneCascadeId

The IDs give us the way to cross the information from the different tables.
The tissue indicate the procedence of the sample, we are only interested on the Blood samples. 
The version_panel indicates the version of the panel used on this sample, so we can know the genes analyzed for each sample. 
The GeneCascadeName is a kind of clasification of the type of disease. 

#### control.samples 
This file contain only the samples belonging to the control group. 

#### all.vars
This file is a list of list and each of this list correspond to one sample, the list contain the information of the variants analyzed for this sample. 

#### all.genes
This file contains the symbol of the genes for each of the panel versions.

The variables it contain are the following, the most important are the bold:

  - id
  - analysisId
  - genomeVersionId
  - **sampleId**
  - **chr**
  - **start**
  - **end**
  - **refBase**
  - **altBase**
  - **uniqueVariantId**
  - depth
  - **freq** indicate the percentage of reads where the variant have been found.
  - **jointQualityScore**

## Data management

We are going to focus only on the samples with tissue equal to blood and DNA (which contain the GCAT control samples), which means that from the 2628 total samples we are going to focus in 2142 belonging to this tissue. 
We also want to take into account another filter for the SNPs which is the **jointQualityScore**, that is a measure of quality, we consider valid data the ones with a lowest value than 4.

So, we search in all.samples, all the 2628 samples the ones from Blood and then cross with all.vars and search the list correpondent to this **sampleId**. With this information of the variables we filter it with **jointQualityScore** > 4, and with a auxiliar function we obtain the genotypes from the **freq** variable. Finally, we 
keep the variables we are interested in: **sampleId**, **chr**, **start**, **end**, **refBase**, **altBase**, **uniqueVariantId** and **genotype.

We repeat this process for all the samples and then make a merge to obtain the data in the following format.

| chromosome |   start  |    end   | reference base | alternative base | unique variant id | sample1 | sample2 | sample3 | ... |
|:----------:|:--------:|:--------:|:--------------:|:----------------:|:-----------------:|:-------:|:-------:|:-------:|:---:|
|      1     | 17****38 | 17****38 |        A       |         G        |        7**5       |   1/1   |   1/1   |   1/1   | ... |
|      1     | 17****88 | 17****88 |        T       |         G        |        7**6       |   1/1   |   1/0   |   1/0   | ... |
|      1     | 17****87 | 17****87 |        T       |         C        |        7**7       |   1/1   |   1/1   |   1/1   | ... |
|     ...    |    ...   |    ...   |       ...      |        ...       |        ...        |   ...   |   ...   |   ...   | ... |

This table have 35893 rows, which is the number of SNPs analyzed, and 2246 columns which are 6 variables + 2431 samples. This data contains 1957 cases and 238 controls.

We have an important variable that is the clasification of the samples, we have the following samples for each clasification (cascade.polides):

| All genes | **Breast** | HBOC | Melanoma | **Ovary** | BC/OC/HBOC + Other | **Control** | HNPCC | Other | Polyposis |
|-----------|:------:|-----:|----------|-------|--------------------|---------|-------|-------|-----------|
| 41        |   813  |  199 | 58       | 292   | 48                 | 238     | 190   | 118   | 198       |

We also create another file which contain the key information of the samples, which will help us to know the clasification and the type of sample.

| analysisId | sampleId | version_panel | individual progeny id  | Gene Cascades Name | Type    | Cascades.polides |
|------------|:--------:|--------------:|------------------------|--------------------|---------|------------------|
| 2**1       |   2**2   |           2.1 | Gi-***4                | ICOCore            | Control | NA               |
| 2**0       |   2**1   |           2.3 | BI-****-16             | MAMA_v1_11-11-16   | Case    | Breast           |
| 2**4       |   2**5   |           2.0 | I-****-05              | RENAL_v1_11-11-16  | Case    | Other            |
| ...        | ...      | ...           | ...                    | ...                | ...     | ...              |

We obtain also from the gene symbol the general gene information and save this information on a table called genes_info, which have the following format. 

| hgnc_symbol | ensembl_gene_id | chromosome_name | start_position | end_position | strand |
|-------------|:---------------:|----------------:|----------------|--------------|--------|
| A2ML1       | ENSG00000166535 |              12 | 8822472        | 8887001      | +      |
| AIP         | ENSG00000110711 |              11 | 67483041       | 67491103     | +      |
| ALK         | ENSG00000171094 |               2 | 29192774       | 29921566     | -      |
| ...         | ...             | ...             | ...            | ...          | ...    |

Then, we also check the genes on the different panels, obtaining a table with the following format:

| genes | panel3 | panel2 | panel1 | panel0 |
|-------|:------:|-------:|--------|--------|
| BAP1  |   2.3  |    2.2 | 2.1    | 2      |
| BARD1 |   2.3  |    2.2 | 2.1    | 2      |
| SDHA  |   2.3  |     NA | NA     | NA     |
| PMS1  | NA     | 2.2    | 2.1    | 2      |

With this we can know if this panel have the gene or not (if it have NA).


## Data visualization 

We make a histogram in order to see the number of genes for each chromosome: 

![Histogram](graphs/genes_per_chr.png)

Finally we use karyoploteR library of R to obtain in a visual way the genes position in each chromosome, this can help us to see if there are some kind of "clusters" of genes.

![Chromosome](graphs/Chr1_genes.png)
![Chromosome](graphs/Chr2_genes.png)
![Chromosome](graphs/Chr3_genes.png)
![Chromosome](graphs/Chr4_genes.png)
![Chromosome](graphs/Chr5_genes.png)
![Chromosome](graphs/Chr6_genes.png)
![Chromosome](graphs/Chr7_genes.png)
![Chromosome](graphs/Chr8_genes.png)
![Chromosome](graphs/Chr9_genes.png)
![Chromosome](graphs/Chr10_genes.png)
![Chromosome](graphs/Chr11_genes.png)
![Chromosome](graphs/Chr12_genes.png)
![Chromosome](graphs/Chr13_genes.png)
![Chromosome](graphs/Chr14_genes.png)
![Chromosome](graphs/Chr15_genes.png)
![Chromosome](graphs/Chr16_genes.png)
![Chromosome](graphs/Chr17_genes.png)
![Chromosome](graphs/Chr18_genes.png)
![Chromosome](graphs/Chr19_genes.png)
![Chromosome](graphs/Chr22_genes.png)
![Chromosome](graphs/ChrX_genes.png)


