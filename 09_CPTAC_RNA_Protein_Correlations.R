library(tidyverse)
library(limma)
library(Biobase)
library(sva)
library(patchwork)

# Reading in the data from the CPTAC data set (Proteomic Data Commons ID: PDC000270)
# Retrieved from LinkedOmics (https://www.linkedomics.org/data_download/CPTAC-PDAC/). 
tumor_rna <- read_tsv("inputs/CPTAC/mRNA_RSEM_UQ_log2_Tumor.tsv", )
tumor_protein <- read_tsv("inputs/CPTAC/proteomics_gene_level_MD_abundance_tumor.tsv", )
tumor_pheno <- read.csv("inputs/CPTAC/molecular_pheno_tumor.csv")

# Turning into data frames
tumor_rna <- as.data.frame(tumor_rna, row.names = 0)
tumor_protein <- as.data.frame(tumor_protein, row.names = 0)

# Renaming columns
colnames(tumor_rna)[1] <- "Genes"
colnames(tumor_protein)[1] <- "Genes"

# Setting row names as genes
rownames(tumor_rna) <- tumor_rna[,1]
rownames(tumor_protein) <- tumor_protein[,1]

# Deleting genes column
tumor_rna <- tumor_rna[,-1]
tumor_protein <- tumor_protein[,-1]


# Making function to get correlation of RNA and protein expression for any given gene
getCor <- function(the_gene = "CEACAM5", cortype = "spearman"){
  gene_rna <- tumor_rna[the_gene,]
  gene_protein <- tumor_protein[the_gene,]
  gene_both <- rbind(gene_rna, gene_protein)
  rownames(gene_both)[1] <- paste0(rownames(gene_rna)[1],"_RNA")
  rownames(gene_both)[2] <- paste0(rownames(gene_protein)[1],"_PROTEIN")
  gene_both <- t(gene_both)
  gene_both <- as.data.frame(gene_both)
  gene_both <- na.omit(gene_both)
  the_cor <- cor(gene_both[,1], gene_both[,2], method = cortype)
  return(the_cor)
}

# Making a list of all the genes that are shared between the RNA and protein data sets
gene_list <- intersect(rownames(tumor_protein), rownames(tumor_rna))

# Getting a list of the spearman correlations between protein-RNA for each gene common to both data sets
cor_list <- list()
for(x in 1:length(gene_list)){
  cor_list[[x]] <- getCor(the_gene = gene_list[x], cortype = "spearman")
  names(cor_list)[[x]] <- gene_list[x]
}

# Ranking the genes by their RNA-protein spearman correlations
cor_vec <- unlist(cor_list)
cor_vec <- na.omit(cor_vec)
sorted_spearman <- sort(cor_vec, decreasing = T)

# Getting a list of the pearson correlations between protein-RNA for each gene common to both data sets
cor_list <- list()
for(x in 1:length(gene_list)){
  cor_list[[x]] <- getCor(the_gene = gene_list[x], cortype = "pearson")
  names(cor_list)[[x]] <- gene_list[x]
}

# Ranking the genes by their RNA-protein pearson correlations
cor_vec <- unlist(cor_list)
cor_vec <- na.omit(cor_vec)
sorted_pearson <- sort(cor_vec, decreasing = T)

# Importing VI-PDAC DEGs
vipdac <- read.csv("outputs/plots_tables_objects/Differential_Expression_VIvsPDAC_all_genes.csv")
rownames(vipdac) <- vipdac$rn
sig_genes <- vipdac  %>% filter((logFC > 0.58 & adj.P.Val < 0.05) | (logFC < -0.58 & adj.P.Val < 0.05))

# Sorting VI-PDAC DEGs by their RNA-Protein correlations
thespearmans <- as.data.frame(sorted_spearman[which(names(sorted_spearman)%in% rownames(sig_genes))])
thepearsons <- as.data.frame(sorted_pearson[which(names(sorted_pearson)%in% rownames(sig_genes))])
colnames(thespearmans)[1] <- "RNA_protein_spearman"
colnames(thepearsons)[1] <- "RNA_protein_pearson"

VI_PDAC <- merge(sig_genes, thespearmans, by = 0, all = T)
rownames(VI_PDAC) <- VI_PDAC$Row.names
VI_PDAC <- VI_PDAC[,-1]
spearman_table_VI_PDAC <- VI_PDAC[order(VI_PDAC$RNA_protein_spearman, decreasing = T),]
VI_PDAC <- merge(spearman_table_VI_PDAC, thepearsons, by = 0, all = T)
pearson_table_VI_PDAC <- VI_PDAC[order(VI_PDAC$RNA_protein_pearson, decreasing = T),]
rownames(pearson_table_VI_PDAC) <- pearson_table_VI_PDAC$Row.names
VI_PDAC <- pearson_table_VI_PDAC[,-1]

write.csv(VI_PDAC, "outputs/plots_tables_objects/spearman_pearson_table_VI_PDAC.csv")
