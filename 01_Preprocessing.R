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
sessionInfo()

# PRE PROCESSING ####

rm(list = ls())

# Read in script containing custom functions
source("scripts/00_Custom_Functions.R")

# This is the ProbeQC filtered data (probeQC settings = geomean probe in all segments / geomean probes within target <= 0.1 & Fails Grubbs outlier test in >= 20% of segments). This was performed in the GeoMx portal.
merged <- read.csv("inputs/ProbeQC_merged_batches.csv", row.names = 1)

# Importing metadata
merged_metadata <- read.csv("inputs/initial_metadata.csv", row.names = 1)

# Renaming all of the tissue types
merged_metadata$Batch <- as.character(merged_metadata$Batch)
merged_metadata$Type[merged_metadata$Type == 'STROMA'] <- 'PDAC'
merged_metadata$Type[merged_metadata$Type == 'NORMAL'] <- 'ND'
merged_metadata$Type[merged_metadata$Type == 'LVI'] <- 'VI'
merged_metadata$Type <- factor(merged_metadata$Type, levels=c('ND','PDAC', 'PNI','VI'))
merged_metadata$Cancer <- "ND"
merged_metadata$Cancer[merged_metadata$Type != 'ND'] <- "COMBINED_PDAC"
merged_metadata$Cancer <- factor(merged_metadata$Cancer, levels=c('COMBINED_PDAC', "ND"))

# Renaming the patients
merged_metadata$Patient.Alias <- NA_character_
merged_metadata$Patient.Alias[merged_metadata$Slide=='2M'] <- "1"
merged_metadata$Patient.Alias[merged_metadata$Slide=='2N'] <- "2"
merged_metadata$Patient.Alias[merged_metadata$Slide=="3K"] <- "3"
merged_metadata$Patient.Alias[merged_metadata$Slide=="4L"] <- "4"
merged_metadata$Patient.Alias[merged_metadata$Slide=="pt238"] <- "5"
merged_metadata$Patient.Alias[merged_metadata$Slide=="pt499"] <- "6"
merged_metadata$Patient.Alias[merged_metadata$Slide=="pt602"] <- "7"
merged_metadata$Patient.Alias[merged_metadata$Slide=="pt63"] <- "8"

#Filtering genes that are above the LOQ in *MORE* than 10 samples. loqfilter is from the custom functions script.
filtered <- loqfilter(merged, 10)
mergedmat <- as.matrix(filtered)
mode(mergedmat) <- "integer"

# saving merged raw data
write.csv(mergedmat, "outputs/preprocessing/mergedmatraw.csv")

## Make dds object
dds <- DESeqDataSetFromMatrix(countData = mergedmat,
                              colData = merged_metadata,
                              design= ~0 + Type)

# using variance stabilizing transformation from DESeq2
vsd <- assay(vst(dds, blind=FALSE))
vsd <- as.data.frame(vsd)

# making heatmap with all of the genes with all of the samples
datScale <- t(apply(vsd,1,scale))
colnames(datScale) <- colnames(vsd)

# FIGURE S1A ---- Heatmap of VST normalized expression data with no batch corrections applied to the data ####
pdf("outputs/preprocessing/pre_corrected_heatmap_with06303.pdf")
Heatmap(datScale, show_row_names = F, show_column_names = T,
        clustering_distance_rows = 'pearson', 
        clustering_distance_columns = 'pearson',
        column_names_gp = grid::gpar(fontsize = 6),
        top_annotation = HeatmapAnnotation(Type =
                                             merged_metadata[,c('Type')],
                                           Patient.Alias = merged_metadata[,c('Patient.Alias')], col = list(Type= c("ND" = "#D6A857", "PDAC" = "#DF2E2E", "PNI" = "#E475FF", "VI" = "#2196F3"),
                                                                                                            Patient.Alias = c("1" = "#00BCD4", "2" = "#000000", "3" = "#CBCBCB", "4" = "#9F5500", "5" = "#7AFF88", "6" = "#3F51B5", "7" = "#FF6600", "8" = "#FFBAE2")),
                                           annotation_legend_param = list(Patient.Alias = list(title = "Patient alias"), Type = list(title = "Tissue type"))))
dev.off()

# Removing the 063_003, due to concern for sample mislabling
data <- merged %>% dplyr::select(-c("Stroma_063_003"))
metadata <- merged_metadata %>% dplyr::slice(-c(which(rownames(merged_metadata) %in%
                                                        c("Stroma_063_003"))))

#Filtering genes that are above the LOQ in *MORE* than 10 samples (with the 063_003 sample removed)
filtered <- loqfilter(data, 10)
mergedmat <- as.matrix(filtered)
mode(mergedmat) <- "integer"

## Making dds object with the 063_003 sample removed
dds <- DESeqDataSetFromMatrix(countData = mergedmat,
                              colData = metadata,
                              design= ~0 + Type)

# using variance stabilizing transformation from DESeq2
vsd <- assay(vst(dds, blind=FALSE))
vsd <- as.data.frame(vsd)

## create the PCA plots for the pre batch corrected data ####

comp_merged<- log2(vsd) %>% prcomp()
table <- subset(comp_merged$rotation, select = c(PC1, PC2))
dtable_merged <- cbind(table, metadata)

# PCA by tissue type
pdf("outputs/preprocessing/pre_corrected_PCA_type_without06303.pdf")
ggplot(dtable_merged, aes_string(x = "PC1", y = "PC2", color = "Type")) +
  geom_point(aes(shape = Type), size = 3) +
  labs(color = "Tissue type", shape = "Tissue type") +
  scale_color_manual(name = "Tissue type", labels = c("ND", "PDAC", "PNI", "VI"), values = c("#D6A857", "#DF2E2E", "#E475FF", "#2196F3"))
dev.off()

# PCA by patient alias
pdf("outputs/preprocessing/pre_corrected_PCA_alias_without06303.pdf")
ggplot(dtable_merged, aes_string(x = "PC1", y = "PC2", color = "Patient.Alias")) +
  geom_point(aes(shape = Type), size = 3) +
  labs(color = "Patient alias", shape = "Tissue type") +
  scale_color_manual(name = "Patient alias", labels = c("1", "2", "3", "4", "5", "6", "7", "8"), values = c("#00BCD4", "#000000", "#CBCBCB", "#9F5500", "#7AFF88", "#3F51B5", "#FF6600", "#FFBAE2"))
dev.off()

## create heatmap of pre batch corrected data ####

# scaling data
datScale <- t(apply(vsd,1,scale))
colnames(datScale) <- colnames(vsd)

pdf("outputs/preprocessing/pre_corrected_heatmap.pdf")
Heatmap(datScale, show_row_names = F, show_column_names = F,
        clustering_distance_rows = 'pearson', 
        clustering_distance_columns = 'pearson',
        top_annotation = HeatmapAnnotation(Type =
                                             metadata[,c('Type')],
                                           Patient.Alias = metadata[,c('Patient.Alias')], col = list(Type= c("ND" = "#D6A857", "PDAC" = "#DF2E2E", "PNI" = "#E475FF", "VI" = "#2196F3"),
                                                                                                     Patient.Alias = c("1" = "#00BCD4", "2" = "#000000", "3" = "#CBCBCB", "4" = "#9F5500", "5" = "#7AFF88", "6" = "#3F51B5", "7" = "#FF6600", "8" = "#FFBAE2")),
                                           annotation_legend_param = list(Patient.Alias = list(title = "Patient alias"), Type = list(title = "Tissue type"))))
dev.off()


# ComBat batch correction ####
# making pheno data and expression set for limma
vsd <- as.matrix(vsd)
pheno = new("AnnotatedDataFrame", data=metadata)
eset = ExpressionSet(assayData = vsd, phenoData = pheno)
edata = Biobase::exprs(eset)

# Applying batch correction for interbatch differences (2 different batches) #####
write.csv(vsd, "outputs/preprocessing/vsd_not_batch_corrected_vsd.csv")

# batch correction on batch
batch_adjusted <- sva::ComBat(dat = edata, batch = metadata$Batch, par.prior = TRUE)

write.csv(batch_adjusted, "outputs/preprocessing/half_batch_corrected_vsd.csv")

# batch correction on patient
batch_adjusted <- sva::ComBat(dat = batch_adjusted, 
                              batch =metadata$Slide, par.prior = TRUE)

# saving data
write.csv(batch_adjusted, "outputs/preprocessing/fully_batch_corrected_vsd.csv")
saveRDS(batch_adjusted, "outputs/preprocessing/batch_adjusted.rds")
write.csv(metadata, "outputs/preprocessing/metadata.csv")
saveRDS(metadata, "outputs/preprocessing/metadata.rds")

comp_merged<- log2(vsd) %>% prcomp()
table <- subset(comp_merged$rotation, select = c(PC1, PC2))
dtable_merged <- cbind(table, metadata)

# PCA plot by tissue type 
pdf("outputs/preprocessing/post_corrected_PCA_type_without06303.pdf")
ggplot(dtable_merged, aes_string(x = "PC1", y = "PC2", color = "Type")) +
  geom_point(aes(shape = Type), size = 3) +
  labs(color = "Tissue type", shape = "Tissue type") +
  scale_color_manual(name = "Tissue type", labels = c("ND", "PDAC", "PNI", "VI"), values = c("#D6A857", "#DF2E2E", "#E475FF", "#2196F3"))
dev.off()

# FIGURE 2B ---- Principal component analysis with patient annotated by color and tissue type annotated by shape #### 
# PCA plot by patient alias
pdf("outputs/preprocessing/post_corrected_PCA_alias_without06303.pdf")
ggplot(dtable_merged, aes_string(x = "PC1", y = "PC2", color = "Patient.Alias")) +
  geom_point(aes(shape = Type), size = 3) +
  labs(color = "Patient alias", shape = "Tissue type") +
  scale_color_manual(name = "Patient alias", labels = c("1", "2", "3", "4", "5", "6", "7", "8"), values = c("#00BCD4", "#000000", "#CBCBCB", "#9F5500", "#7AFF88", "#3F51B5", "#FF6600", "#FFBAE2"))
dev.off()
#####

# MAKE HEATMAP OF BATCH CORRECTED DATA ####

# scaling data + transposing data + renaming columns
datScale <- t(apply(batch_adjusted,1,scale))
colnames(datScale) <- colnames(batch_adjusted)

# heatmap of batch corrected data
pdf("outputs/preprocessing/post_corrected_heatmap.pdf")
Heatmap(datScale, show_row_names = F, show_column_names = T,
        clustering_distance_rows = 'pearson', 
        clustering_distance_columns = 'pearson',
        top_annotation = HeatmapAnnotation(Type =
                                             metadata[,c('Type')],
                                           Patient.Alias = metadata[,c('Patient.Alias')], col = list(Type= c("ND" = "#D6A857", "PDAC" = "#DF2E2E", "PNI" = "#E475FF", "VI" = "#2196F3"),
                                                                                                     Patient.Alias = c("1" = "#00BCD4", "2" = "#000000", "3" = "#CBCBCB", "4" = "#9F5500", "5" = "#7AFF88", "6" = "#3F51B5", "7" = "#FF6600", "8" = "#FFBAE2")),
                                           annotation_legend_param = list(Patient.Alias = list(title = "Patient alias"), Type = list(title = "Tissue type"))),
        column_names_gp = gpar(fontsize = 5))
dev.off()
