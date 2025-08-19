library(filesstrings)
library(edgeR)
library(projectR)
library(purrr)
library(CoGAPS)
library(psych)
library(FSA)
library(tidyverse)
library(sva)
library(limma)
library(Biobase)
library(patchwork)
library(data.table)
library(ggrepel)
library(DESeq2)
library(gridExtra)
library(caroline)
library('DT')
library('EnhancedVolcano')
library('ComplexHeatmap')
library(webshot)
library(VennDiagram)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(projectR)
library(Matrix)
sessionInfo()

# Reading in custom functions
source("scripts/00_Custom_Functions.R")

# reading in the files for the gallinger bulk RNA seq data
datadir <- "inputs/PanCuRx_rnacount/"
TxtFiles <- dir(datadir, pattern=".txt")
list_of_txtfiles <- list()
# storing each of those text files as a list
list_of_txtfiles <- list()
for(x in 1:length(TxtFiles)){
  print(x)
  list_of_txtfiles[[x]] <- read.csv(paste0(datadir,TxtFiles[x]), sep = "\t", header = F)
}

# joining all of the text files into a single data frame
thedata <- list_of_txtfiles %>% purrr::reduce(full_join, by = "V1") %>% as.data.frame()
colnames(thedata)[2:ncol(thedata)] <- TxtFiles
rownames(thedata) <- thedata[,1]
thedata <- thedata[,-1]

# editing the column names
colnames(thedata) <- gsub("_htseq_count_all.txt", "", colnames(thedata))
thedata <- as.data.frame(thedata)

# filtering for the data present in paired_unpaired (from the paper)
paired_unpaired <- read.csv("inputs/paired_unpaired.csv", header = T)

### this code turns ensembl_gene_ids into hgnc symbols
# ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
# symbols <- getBM(attributes=c("ensembl_gene_id",'hgnc_symbol'), 
#                  filters = 'ensembl_gene_id', 
#                  values = rownames(thedata), 
#                  mart = ensembl)
#saveRDS(ensembl, "inputs/ensembl.rds")
#saveRDS(symbols, "inputs/symbols.rds")
# reading in the ensembl and hgnc symbols
ensembl <- readRDS("inputs/ensembl.rds")
symbols <- readRDS("inputs/symbols.rds")

# filtering data by hgnc symbols, removing duplicate hgnc symbols, removing IDs without corresponding hgnc symbols
thedata$ensembl_gene_id <- rownames(thedata)
newdata <- merge(thedata, symbols, by = "ensembl_gene_id")
newdata <- newdata[-which(duplicated(newdata$hgnc_symbol)),]
newdata <- newdata[-which(newdata$hgnc_symbol == "")]
newdata <- newdata[,-1]
rownames(newdata) <- newdata$hgnc_symbol
newdata <- newdata[, -ncol(newdata)]

# finding out which samples are present in the ICGC data but not in the data that was shared over email

exp_seq <- read.csv("inputs/exp_seq.tsv", sep = "\t")
exp_seq_wts = exp_seq %>% filter(analysis_id == 'PACA_CA-WTS')
exp_seq_wts$gene_id = as.factor(exp_seq_wts$gene_id)
exp_seq_wts$submitted_sample_id = as.factor(exp_seq_wts$submitted_sample_id)

# create matrix
expMat <- as.matrix(sparseMatrix(
  i = as.numeric(exp_seq_wts$gene_id),
  j = as.numeric(exp_seq_wts$submitted_sample_id),
  x = exp_seq_wts$raw_read_count,
  dimnames = list(levels(exp_seq_wts$gene_id),
                  levels(exp_seq_wts$submitted_sample_id))))

ega_idents <- colnames(newdata)
icgc_idents <- colnames(expMat)

left_overs <- icgc_idents[which(!icgc_idents %in% ega_idents)]
left_overs <- expMat %>% as.data.frame() %>%dplyr::select(all_of(left_overs))

write.csv(left_overs,"outputs/left_over_icgc_samples.csv")

