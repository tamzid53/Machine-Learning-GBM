
setwd("/GBM/")

####Related to Data Import, Cleaning and Preprocessing Steps 1-10


##Load R packages:

library(UCSCXenaTools)
library(data.table)
library(R.utils)
library(dplyr)


data(XenaData)
write.csv(XenaData, "00_tblXenaHubInfo.csv")


####[Main Text: Step 5] Select then download target data sets from Xena Data Hubs.
##[Step 5-a] Target=RSEM expected counts provided by the UCSC toil Recompute Compendium

getOption('timeout')
options(timeout=10000)

GeneExpectedCnt_toil =XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TcgaTargetGtex_gene_expected_count")
XenaQuery(GeneExpectedCnt_toil) %>%
  XenaDownload(destdir ="./")   ###, max_try = 1L)




##[Step 5-b] Target =TCGA Clinical data.

paraCohort = "TCGA Glioblastoma"
paraDatasets ="TCGA.GBM.sampleMap/GBM_clinicalMatrix"


Clin_TCGA = XenaGenerate(subset = XenaHostNames =="tcgaHub") %>%
  XenaFilter(filterCohorts = paraCohort) %>%
  XenaFilter(filterDatasets = paraDatasets)
XenaQuery(Clin_TCGA) %>%
  XenaDownload(destdir = "./")


##[Step 5-c] Target = TCGA Survival Data.

Surv_TCGA =XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TCGA_survival_data")
XenaQuery(Surv_TCGA) %>%
  XenaDownload(destdir ="./")



##[Step 5-d] Target = GTEx Phenotype data.

Pheno_GTEx = XenaGenerate(subset =XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TcgaTargetGTEX_phenotype")
XenaQuery(Pheno_GTEx) %>%
  XenaDownload(destdir ="./")




#### [Main Text: Step 6] Subset expression data to include only desired samples.

## [Step 6-a] Retrive IDs for GTEx normal samples of desired tissue type(s).

filterGTEx01 = fread("TcgaTargetGTEX_phenotype.txt.gz")
names(filterGTEx01) = gsub("\\_", "", names(filterGTEx01))



paraStudy = "GTEX" # Setting "GTEx" as the study of interest
paraPrimarySiteGTEx = "Brain" # Setting "brain" as the primary site of interest

filterGTEx02 = subset(filterGTEx01,	study == paraStudy &  primarysite == paraPrimarySiteGTEx )

dim(filterGTEx02)
#[1] 1152    7


## [Step 6-b] Retrive IDs for TCGA primary tumor samples of desired histological type(s).

filterTCGA01 = fread(paraDatasets)
names(filterTCGA01) = gsub("\\_", "", names(filterTCGA01))



paraSampleType = "Primary Tumor" #Setting "Primary Tumor" as the sample type of interest.
paraPrimarySiteTCGA = "Brain" #Setting "Glioblastoma" as the primary site of interest.


filterTCGA02 = subset(filterTCGA01,	sampletype == paraSampleType &	primarysite == paraPrimarySiteTCGA) 
#dim(filterTCGA02)
#[1] 602 129




filterExpr = c(filterGTEx02$sample, filterTCGA02$sampleID, "sample")

ExprSubsetBySamp =fread("TcgaTargetGtex_gene_expected_count.gz", select=filterExpr)





