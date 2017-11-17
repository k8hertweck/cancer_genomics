#### Downloading, inspecting, and parsing data from TCGA ####
# BRCA is breast cancer data

## install TCGA bioconductor tools
# identify location for Bioconductor tools
#source("https://bioconductor.org/biocLite.R")
# install Bioconductor tools
#biocLite("TCGAbiolinks")
#biocLite("SummarizedExperiment")
#biocLite("maftools")
# load Bioconductor tools
library(TCGAbiolinks)
library(SummarizedExperiment)
library(maftools)

# install packages from CRAN
#install.packages("dplyr")
# load packages from CRAN
library(dplyr)

#### Set up project ####

# create directory structure
dir.create("data")
dir.create("figures")

#### Identify TCGA data available #### 

# show all available projects
getGDCprojects()$project_id

# view data available for breast cancer (TCGA-BRCA)
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")

#### Clinical data ####

# download and read clinical data for breast cancer into R
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
# write data to file
write.table(clinical, "data/clinicalBRCA.csv")

# inspecting variables of interest
str(clinical) # 1097 total records
table(clinical$race)
table(clinical$vital_status) 
table(clinical$morphology) 
clinical$days_to_death
clinical$bcr_patient_barcode # patient ID

#### Gene expression (transcriptome) data ####

# identify desired data
query_fpkm <- GDCquery(project = "TCGA-BRCA", 
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification", 
                       workflow.type = "HTSeq - FPKM-UQ")
# download data
# downloads query data into GDCdata/
# data files are rather large; this step takes awhile!
GDCdownload(query_fpkm)
# read downloaded data into R
fpkm <- GDCprepare(query_fpkm)
# the commands above may create the following supplemental files:
# Human_genes_GRCh38_p10.rda
# MANIFEST.txt
# These files can be moved to GDCdata/

# save imported object to file
save(fpkm, file="GDCdata/geneExpressionBRCA.RData")

# load saved data
load("GDCdata/geneExpressionBRCA.RData")

# inspect structure and features of object
assayNames(fpkm)
head(assay(fpkm), 1)
colSums(assay(fpkm))
rowRanges(fpkm)
colData(fpkm) # metadata
colnames(colData(fpkm)) # just metadata column names
# show gene names
rowRanges(fpkm)$external_gene_name

# extract genes of interest
# create object of ensembl_gene_id and external_gene_name
genes <-rowData(fpkm)
# find BRCA1 and BRCA2
brca1 <- grep("brca", genes$external_gene_name, ignore.case = TRUE)
genes[brca1, ]
# BRCA1 ENSG00000012048
# BRCA2 ENSG00000139618

##  assemble dataset for genes of interest and metadata
fpkmDat <- as.data.frame(t(assays(fpkm)[[1]])) # extract expression data
colnames(fpkmDat) # print gene names
rownames(fpkmDat) # show sample names
# extract gene data for target genes
fpkmGene <- fpkmDat %>%
  select(ENSG00000012048, ENSG00000139618)
# extract metadata
metaDat <-as.data.frame(colData(fpkm))
# bind metadata to gene expression data
fpkmGene <- cbind(fpkmGene, metaDat)
# create object of gene names in order
geneNames <- c("BRCA1", "BRCA2")
# create object of metadata names
metaNames <- colnames(colData(fpkm))
# replace column names
colnames(fpkmGene) <- c(geneNames, metaNames)

## clean data
# remove troublesome metadata
fpkmGene <- select(fpkmGene, -treatments)
# save aggregated data to file
write.table(fpkmGene, "data/targetGeneBRCA.csv")

#### Simple Nucelotide Variation ####

#https://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html

# download somatic mutations as maf files (save as csv)
query_maf_muse <- GDCquery_Maf("BRCA", save.csv =TRUE, pipelines = "muse")
query_maf_varscan2 <- GDCquery_Maf("BRCA", save.csv =TRUE, pipelines = "varscan2")
query_maf_somaticsniper <- GDCquery_Maf("BRCA", save.csv =TRUE, pipelines = "somaticsniper")
query_maf_mutect2 <- GDCquery_Maf("BRCA", save.csv =TRUE, pipelines = "mutect2")

# create maf object (without clinical data)
maf_muse <- read.maf(query_maf_muse, useAll = FALSE)
# visual summary of data
plotmafSummary(maf = maf_muse, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
# draw oncoplot
oncoplot(maf = maf_muse, top = 10, removeNonMutated = TRUE)
# assess transitions and transversions
maf_titv <- titv(maf = maf_muse, plot = FALSE, useSyn = TRUE)
# plot transitions and transversions
plotTiTv(res = maf_titv)

# extract barcodes from maf query
barcodes <- sort(query_maf$Tumor_Sample_Barcode) %>%
  unique
str_trunc(barcodes, 12, side = "right", ellipsis = "")

# extract barcodes from clinical data
sub_id <- sort(clinical$submitter_id) %>%
  unique()
clinical <- rename(clinical, Tumor_Sample_Barcode = submitter_id)

# download all snp data (including germline)
query_snp <- GDCquery(project = "TCGA-BRCA", 
                      data.category = "Simple Nucleotide Variation")


# create maf object with clinical data attached
maf <- read.maf(query_maf, clinicalData = clinical, useAll = FALSE)

