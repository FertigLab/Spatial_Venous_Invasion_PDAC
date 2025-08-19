library(CoGAPS)
library(GSVA)
library(tidyverse)
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
library(msigdbr)
library(fgsea)
library(GSVA)
sessionInfo()

source("scripts/00_Custom_Functions.R")

# reading in the VI subtype annotations
annots <- read.csv("inputs/annotations.csv", header= F)
colnames(annots) <- c("Slide", "Number", "Classification")
annots <- as.data.frame(annots)
thedata <- read.csv("outputs/preprocessing/fully_batch_corrected_vsd.csv", row.names = 1)
metadata <- read.csv("outputs/preprocessing/metadata.csv", row.names = 1)
# removing stroma 063_003
annots <- annots[-87,]
metadata <- cbind(metadata, annots$Classification)
colnames(metadata)[16] <- "VI_Classification"
rownames(metadata) <- metadata$Name


# Volcano plots of each of the subtypes
# reading in the fully batch corrected data
thedata <- read.csv("outputs/preprocessing/fully_batch_corrected_vsd.csv", row.names = 1)

metadata$mixed_classification <- "NA"
metadata$mixed_classification[which(metadata$Type == "PDAC")] <- "PDAC"
metadata$mixed_classification[which(metadata$VI_Classification == "A")] <- "A"
metadata$mixed_classification[which(metadata$VI_Classification == "B")] <- "B"
metadata$mixed_classification[which(metadata$VI_Classification == "C")] <- "C"

thedata <- thedata[,metadata$mixed_classification != "NA"]
metadata <- metadata[metadata$mixed_classification != "NA",]

thedata <- thedata[,metadata$mixed_classification != "C"]
metadata <- metadata[metadata$mixed_classification != "C",]

#theviindex <- which(annots$Classification != "X")
#thedata <- thedata[,theviindex]
#metadata <- metadata[theviindex,]
pheno = new("AnnotatedDataFrame", data=metadata)

# Making matrix design
design = model.matrix(~0 + pheno$mixed_classification)
colnames(design) <- c("A", "B", "PDAC")

# Making contrasts for differential expression
contr.matrix <- makeContrasts(
  BvsA = B - A,
  BvsPDAC = B-PDAC,
  AvsPDAC = A-PDAC,
  levels = colnames(design))

# doing limma differential expression on the batch adjusted matrix
eset = ExpressionSet(assayData = as.matrix(thedata), phenoData = pheno)
fit <- lmFit(eset, design)
vfit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(vfit)
tfit <-treat(efit)
results_efit<-decideTests(efit, adj.P.val=0.05)
summary_efit <- summary(results_efit)

# performing the individual comparisons
B.vs.A<- topTreat(efit, coef=1, n=Inf)
B.vs.PDAC<- topTreat(efit, coef=2, n=Inf)
A.vs.PDAC<- topTreat(efit, coef=3, n=Inf)

# making a table of the p values for all the comparisons
inlike_destructive <- data.frame("Genes" = rownames(B.vs.A), "group_1" = "VI-IN-Like", "group_2" = "VI-Destructive", "p" = B.vs.A$adj.P.Val, "LFC" = B.vs.A$logFC)
inlike_pdac <- data.frame("Genes" = rownames(B.vs.PDAC), "group_1" = "VI-IN-Like", "group_2" = "PDAC", "p" = B.vs.PDAC$adj.P.Val, "LFC" = B.vs.PDAC$logFC)
destructive_pdac <- data.frame("Genes" = rownames(A.vs.PDAC), "group_1" = "VI-Destructive", "group_2" = "PDAC", "p" = A.vs.PDAC$adj.P.Val, "LFC" = A.vs.PDAC$logFC)
p_table_2 <- rbind(inlike_destructive, inlike_pdac, destructive_pdac)
saveRDS(p_table_2, "outputs/plots_tables_objects/p_table_2.rds")

# FIGURE S8A ---- making volcano plots for VI-IN-Like vs VI-Destructive ####
# making volcano plots for each comparisons
volcanoPlot(B.vs.A, title_1 = "VI-IN-Like", title_2 = "VI-Destructive", upcolor = "#36DAFF", downcolor = "#9C27B0", batch = "vsd", directory = "outputs/plots_tables_objects/")
#####

# reading in the cogaps pattern markers
jp4thresh <- readRDS('outputs/thresholded_patterns.rds')

# getting the hallmark gene sets
all_gene_sets = msigdbr(species = "human", category = "H")
msigdbr_list = split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name)

# performing fgsea with logFC rank
the_rank <- B.vs.A$logFC
names(the_rank) <- rownames(B.vs.A)
BvsAfgsea <- fgsea(pathways = msigdbr_list, stats = the_rank)

# making lollipop plot with fgsea output from logFC rank
BvsAfgsea <- as.data.frame(BvsAfgsea)
NES <- BvsAfgsea[order(BvsAfgsea$NES, decreasing = T),]

# making lollipop plot
positives <- NES %>% filter(NES > 0)
negatives <- NES %>% filter(NES < 0)

# making lollipop plot
positives <- positives[order(positives$NES, decreasing = F), ]
negatives <- negatives[order(negatives$NES, decreasing = F), ]

# makign lollipop plot
positives <- positives %>% mutate(group = "positives")
negatives <- negatives %>% mutate(group = "negatives")

# making lollipop plot
thedata <- rbind(negatives, positives)
theleadingedge <- t(plyr::ldply(thedata$leadingEdge, rbind))
colnames(theleadingedge) <- thedata$pathway[!is.na(thedata$pathway)]
thedata <- thedata %>% dplyr::select("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "group")
thedata$pathway <- factor(thedata$pathway, levels = thedata$pathway)
thedatatrimmed <- thedata %>% filter(NES >= 2 | NES <= -2)

write.csv(thedata, "outputs/plots_tables_objects/BvsA_fgsea_results.csv")
write.csv(theleadingedge, "outputs/plots_tables_objects/BvsA_fgsea_leadingedge.csv")

# FIGURE 6B ---- making lollipop plot of GSEA for VI-IN-Like vs Destructive ####
pdf('outputs/plots_tables_objects/IN-LikevsDestructive_fgsea_LFC_lollipop.pdf')
ggplot(thedata, aes(x=pathway, y=NES)) +
  geom_segment(aes(x=pathway, xend=pathway, y=0, yend=NES), size=0.5, alpha=0.5) +
  geom_point(size = 3, aes(fill= group), alpha=1, shape=21, stroke=1.4) +
  scale_fill_manual(values = c("#9C27B0", "#36DAFF")) +
  theme_light() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
  ) +
  xlab("") +
  ylab("Normalized Enrichment Score") +
  coord_flip() +
  theme(axis.text.y = element_text(size = 6.8)) +
  scale_colour_manual("red", "blue")
dev.off()
#####

# FIGURE S6A ---- barcode plot for EMT gene sets ####
pdf("outputs/plots_tables_objects/emt_LFC_barcodeplot.pdf")
index_emt <- as.vector(na.omit(match(msigdbr_list$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, rownames(B.vs.A))))
logfc <- B.vs.A$logFC
# setting quantiles to be extremely large so that the barcode plot isn't arbitrarily colored blue/grey/red. 
barcodeplot(logfc, index_emt, labels = c("VI-Destructive", "VI-IN-Like"), quantiles = c(-9999,9999)*sqrt(9999))
dev.off()
#####

# FIGURE S6B ---- barcode plot for MYC gene set ####
pdf("outputs/plots_tables_objects/myc_LFC_barcodeplot.pdf")
index_mycv1 <- as.vector(na.omit(match(msigdbr_list$HALLMARK_MYC_TARGETS_V1, rownames(B.vs.A))))
logfc <- B.vs.A$logFC
# setting quantiles to be extremely large so that the barcode plot isn't arbitrarily colored blue/grey/red. 
barcodeplot(logfc, index_mycv1, labels = c("VI-Destructive", "VI-IN-Like"), quantiles = c(-9999,9999)*sqrt(9999))
dev.off()
#####

# getting barcode plots for Classical and Basal
classical <- c('BTNL8','FAM3D','ATAD4','AGR3','CTSE','LOC400573','LYZ','TFF2',
               'TFF1','ANXA10','LGALS4','PLA2G10','CEACAM6','VSIG2','TSPAN8',
               'ST6GALNAC1','AGR2','TFF3','CYP3A7','MYO1A','CLRN3','KRT20',
               'CDH17','SPINK4','REG4')

basal <- c('VGLL','UCA1','S100A2','LY6D','SPRR3','SPRR1B','LEMD1','KRT15',
           'CTSL2','DHRS9','AREG','CST6','SERPINB3','KRT6C','KRT6A','SERPINB4',
           'FAM83A','SCEL','FGFBP1','KRT7','KRT17','GPR87','TNS4','SLC2A1',
           'ANXA8L2')

# FIGURE 6C - barcode plots for classical/basal signatures like in Destructive vs IN-Like ####
pdf("outputs/plots_tables_objects/classical_LFC_barcodeplot.pdf")
index_classical <- as.vector(na.omit(match(classical, rownames(B.vs.A))))
logfc <- B.vs.A$logFC
# setting quantiles to be extremely large so that the barcode plot isn't arbitrarily colored blue/grey/red. 
barcodeplot(logfc, index_classical, labels = c("VI-Destructive", "VI-IN-Like"), quantiles = c(-9999,9999)*sqrt(9999))
dev.off()

pdf("outputs/plots_tables_objects/basal_LFC_barcodeplot.pdf")
index_basal <- as.vector(na.omit(match(basal, rownames(B.vs.A))))
logfc <- B.vs.A$logFC
# setting quantiles to be extremely large so that the barcode plot isn't arbitrarily colored blue/grey/red. 
barcodeplot(logfc, index_basal, labels = c("VI-Destructive", "VI-IN-Like"), quantiles = c(-9999,9999)*sqrt(9999))
dev.off()
#####

# getting barcode plots for mesenchymal and collective in Destructive vs IN-Like
collective <- c("TM4SF4", "ONECUT2", "RND1", "GATA6", "HNF1A")

mesenchymal <- c("ADAMTS9", "PTGS2", "TGFBI", "SERPINH1", "SPARC", "MFAP2", "MXRA8")

# FIGURE S6C - barcode plots for collective/mesenchymal in Destructive vs IN-Like ####
pdf("outputs/plots_tables_objects/collective_LFC_barcodeplot.pdf")
index_collective <- as.vector(na.omit(match(collective, rownames(B.vs.A))))
logfc <- B.vs.A$logFC
# setting quantiles to be extremely large so that the barcode plot isn't arbitrarily colored blue/grey/red. 
barcodeplot(logfc, index_collective, labels = c("VI-Destructive", "VI-IN-Like"), quantiles = c(-9999,9999)*sqrt(9999))
dev.off()

pdf("outputs/plots_tables_objects/mesenchymal_LFC_barcodeplot.pdf")
index_mesenchymal <- as.vector(na.omit(match(mesenchymal, rownames(B.vs.A))))
logfc <- B.vs.A$logFC
# setting quantiles to be extremely large so that the barcode plot isn't arbitrarily colored blue/grey/red. 
barcodeplot(logfc, index_mesenchymal, labels = c("VI-Destructive", "VI-IN-Like"), quantiles = c(-9999,9999)*sqrt(9999))
dev.off()
#####

# CoGAPS pattern plots for Destructive vs IN-Like
 metadata <- read.csv("outputs/preprocessing/metadata.csv", row.names = 1)

# removing stroma 063_003 from the annotation annots <- annots[-87,]
metadata <- cbind(metadata, annots$Classification)
colnames(metadata)[16] <- "VI_Classification"
rownames(metadata) <- metadata$Name


# reading in the cogaps object
cogaps <- readRDS("inputs/batch_merged_cogapsP4.rds")

# renaming the different tissue types
newmetadata <- metadata
newmetadata$VI_Classification[which(newmetadata$VI_Classification == "X")] <- "ND"
newmetadata$VI_Classification[which(newmetadata$VI_Classification == "A")] <- "VI-Destructive"
newmetadata$VI_Classification[which(newmetadata$VI_Classification == "B")] <- "VI-IN-Like"
newmetadata$VI_Classification[which(newmetadata$VI_Classification == "C")] <- "VI-Conventional"
newmetadata$VI_Classification[which(newmetadata$Type == "PDAC")] <- "PDAC"
newmetadata$VI_Classification[which(newmetadata$Type == "NORMAL")] <- "ND"
newmetadata$VI_Classification[which(newmetadata$Type == "PNI")] <- "PNI"
saveRDS(newmetadata, "outputs/preprocessing/metadata_with_VI_annotations.rds" )
write.csv(newmetadata, "outputs/preprocessing/metadata_with_VI_annotations.csv" )
newmetadata <- cbind(newmetadata,cogaps@sampleFactors)

# filtering in and factoring the tissue types of interest
newmetadata <- newmetadata %>% filter(VI_Classification == "ND" | VI_Classification == "VI-Destructive" | VI_Classification == "VI-IN-Like" | VI_Classification == "PDAC")
newmetadata$VI_Classification <- factor(newmetadata$VI_Classification, levels = c("ND", "VI-IN-Like", "VI-Destructive", "PDAC"))

# isolating the pattern scores for each tissue type.. this is for getting stats
pdacpat1 <- newmetadata %>% filter(VI_Classification == "PDAC") %>% dplyr::select("Pattern_1") %>% as.vector() %>% unlist() %>% as.vector()
vdpat1 <- newmetadata %>% filter(VI_Classification == "VI-Destructive") %>% dplyr::select("Pattern_1") %>% as.vector() %>% unlist() %>% as.vector()
vinpat1 <- newmetadata %>% filter(VI_Classification == "VI-IN-Like") %>% dplyr::select("Pattern_1") %>% as.vector() %>% unlist() %>% as.vector()

pdacpat2 <- newmetadata %>% filter(VI_Classification == "PDAC") %>% dplyr::select("Pattern_2") %>% as.vector() %>% unlist() %>% as.vector()
vdpat2 <- newmetadata %>% filter(VI_Classification == "VI-Destructive") %>% dplyr::select("Pattern_2") %>% as.vector() %>% unlist() %>% as.vector()
vinpat2 <- newmetadata %>% filter(VI_Classification == "VI-IN-Like") %>% dplyr::select("Pattern_2") %>% as.vector() %>% unlist() %>% as.vector()

pdacpat3 <- newmetadata %>% filter(VI_Classification == "PDAC") %>% dplyr::select("Pattern_3") %>% as.vector() %>% unlist() %>% as.vector()
vdpat3 <- newmetadata %>% filter(VI_Classification == "VI-Destructive") %>% dplyr::select("Pattern_3") %>% as.vector() %>% unlist() %>% as.vector()
vinpat3 <- newmetadata %>% filter(VI_Classification == "VI-IN-Like") %>% dplyr::select("Pattern_3") %>% as.vector() %>% unlist() %>% as.vector()

pdacpat4 <- newmetadata %>% filter(VI_Classification == "PDAC") %>% dplyr::select("Pattern_4") %>% as.vector() %>% unlist() %>% as.vector()
vdpat4 <- newmetadata %>% filter(VI_Classification == "VI-Destructive") %>% dplyr::select("Pattern_4") %>% as.vector() %>% unlist() %>% as.vector()
vinpat4 <- newmetadata %>% filter(VI_Classification == "VI-IN-Like") %>% dplyr::select("Pattern_4") %>% as.vector() %>% unlist() %>% as.vector()

# doing stats between tissue comparisons
pdac_vd1 <- wilcox.test(pdacpat1, vdpat1)
pdac_vin1 <- wilcox.test(pdacpat1, vinpat1)
vin_vd1 <- wilcox.test(vinpat1, vdpat1)

pdac_vd2 <- wilcox.test(pdacpat2, vdpat2)
pdac_vin2 <- wilcox.test(pdacpat2, vinpat2)
vin_vd2 <- wilcox.test(vinpat2, vdpat2)

pdac_vd3 <- wilcox.test(pdacpat3, vdpat3)
pdac_vin3 <- wilcox.test(pdacpat3, vinpat3)
vin_vd3 <- wilcox.test(vinpat3, vdpat3)

pdac_vd4 <- wilcox.test(pdacpat4, vdpat4)
pdac_vin4 <- wilcox.test(pdacpat4, vinpat4)
vin_vd4 <- wilcox.test(vinpat4, vdpat4)

# making dataframes of all the p values for each tissue type and comparison
pdac_vd1p <- data.frame("Pattern" = "Pattern_1", "group_1" = "PDAC", "group_2" = "VI-Destructive", "p" = pdac_vd1$p.value)
pdac_vin1p <- data.frame("Pattern" = "Pattern_1", "group_1" = "PDAC", "group_2" = "VI-IN-Like", "p" = pdac_vin1$p.value)
vin_vd1p <- data.frame("Pattern" = "Pattern_1", "group_1" = "VI-IN-Like", "group_2" = "VI-Destructive", "p" = vin_vd1$p.value)

pdac_vd2p <- data.frame("Pattern" = "Pattern_2", "group_1" = "PDAC", "group_2" = "VI-Destructive", "p" = pdac_vd2$p.value)
pdac_vin2p <- data.frame("Pattern" = "Pattern_2", "group_1" = "PDAC", "group_2" = "VI-IN-Like", "p" = pdac_vin2$p.value)
vin_vd2p <- data.frame("Pattern" = "Pattern_2", "group_1" = "VI-IN-Like", "group_2" = "VI-Destructive", "p" = vin_vd2$p.value)

pdac_vd3p <- data.frame("Pattern" = "Pattern_3", "group_1" = "PDAC", "group_2" = "VI-Destructive", "p" = pdac_vd3$p.value)
pdac_vin3p <- data.frame("Pattern" = "Pattern_3", "group_1" = "PDAC", "group_2" = "VI-IN-Like", "p" = pdac_vin3$p.value)
vin_vd3p <- data.frame("Pattern" = "Pattern_3", "group_1" = "VI-IN-Like", "group_2" = "VI-Destructive", "p" = vin_vd3$p.value)

pdac_vd4p <- data.frame("Pattern" = "Pattern_4", "group_1" = "PDAC", "group_2" = "VI-Destructive", "p" = pdac_vd4$p.value)
pdac_vin4p <- data.frame("Pattern" = "Pattern_4", "group_1" = "PDAC", "group_2" = "VI-IN-Like", "p" = pdac_vin4$p.value)
vin_vd4p <- data.frame("Pattern" = "Pattern_4", "group_1" = "VI-IN-Like", "group_2" = "VI-Destructive", "p" = vin_vd4$p.value)

# making table of p values
p_table <- rbind(pdac_vd1p, pdac_vin1p, vin_vd1p, pdac_vd2p, pdac_vin2p, vin_vd2p, pdac_vd3p, pdac_vin3p, vin_vd3p, pdac_vd4p, pdac_vin4p, vin_vd4p)

# FIGURE 6D ---- boxplot with CoGAPS patterns across the tissue types ####
newmetadatapat <- newmetadata %>% dplyr::select("VI_Classification", "Pattern_1")
data_for_p_val_manual <- p_table %>% filter(Pattern == "Pattern_1") %>% dplyr::select("group_1", "group_2", "p")
#data_for_p_val_manual$label <- paste0("p = ", signif(data_for_p_val_manual$p, digits = 3))
#data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "label")
colnames(data_for_p_val_manual) <- c("group1", "group2", "p")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)

# Convert p-values to asterisks
datman$p.signif <- sapply(datman$p, p_to_asterisk)

#datman$p <- signif(datman$p,digits=3)
yposition <- max(newmetadatapat[,"Pattern_1"]) * 1.05
stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
pdf("outputs/plots_tables_objects/Pattern_1_boxplot.pdf", width = 4, height = 6)
ggplot(newmetadata,aes(x = VI_Classification, y = Pattern_1, color = VI_Classification)) +
  geom_boxplot() + 
  geom_jitter() + 
  labs(color = "", x = "", y = "Pattern weight", title = "Pattern 1") +
  stat.test +
  theme(axis.text.x = element_text(color = "black",
                                   size = 12, angle = 45, hjust = 1),
        plot.title = element_text(color = "black", size = 18, hjust = 0.5)) +
  scale_color_manual(name = "Tissue type", labels = c("ND", "VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#D6A857", "#36DAFF", "#9C27B0", "#DF2E2E"))
dev.off()

# making boxplot with pattern 2 across the tissue types, with statistics 
newmetadatapat <- newmetadata %>% dplyr::select("VI_Classification", "Pattern_2")
data_for_p_val_manual <- p_table %>% filter(Pattern == "Pattern_2") %>% dplyr::select("group_1", "group_2", "p")
# data_for_p_val_manual$label <- paste0("p = ", signif(data_for_p_val_manual$p, digits = 3))
# data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "label")
colnames(data_for_p_val_manual) <- c("group1", "group2", "p")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)

# Convert p-values to asterisks
datman$p.signif <- sapply(datman$p, p_to_asterisk)

#datman$p <- signif(datman$p,digits=3)
yposition <- max(newmetadatapat[,"Pattern_2"]) * 1.05
stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
pdf("outputs/plots_tables_objects/Pattern_2_boxplot.pdf", width = 4, height = 6)
ggplot(newmetadata,aes(x = VI_Classification, y = Pattern_2, color = VI_Classification)) +
  geom_boxplot() + 
  geom_jitter() + 
  labs(color = "", x = "", y = "Pattern weight", title = "Pattern 2") +
  stat.test +
  theme(axis.text.x = element_text(color = "black",
                                   size = 12, angle = 45, hjust = 1),
        plot.title = element_text(color = "black", size = 18, hjust = 0.5)) +
  scale_color_manual(name = "Tissue type", labels = c("ND", "VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#D6A857", "#36DAFF", "#9C27B0", "#DF2E2E"))

dev.off()

# making boxplot with pattern 3 across the tissue types, with statistics 

newmetadatapat <- newmetadata %>% dplyr::select("VI_Classification", "Pattern_3")
data_for_p_val_manual <- p_table %>% filter(Pattern == "Pattern_3") %>% dplyr::select("group_1", "group_2", "p")
# data_for_p_val_manual$label <- paste0("p = ", signif(data_for_p_val_manual$p, digits = 3))
# data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "label")
colnames(data_for_p_val_manual) <- c("group1", "group2", "p")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)

# Convert p-values to asterisks
datman$p.signif <- sapply(datman$p, p_to_asterisk)

#datman$p <- signif(datman$p,digits=3)
yposition <- max(newmetadatapat[,"Pattern_3"]) * 1.05
stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
pdf("outputs/plots_tables_objects/Pattern_3_boxplot.pdf", width = 4, height = 6)
ggplot(newmetadata,aes(x = VI_Classification, y = Pattern_3, color = VI_Classification)) +
  geom_boxplot() + 
  geom_jitter() + 
  labs(color = "", x = "", y = "Pattern weight", title = "Pattern 3") +
  stat.test +
  theme(axis.text.x = element_text(color = "black",
                                   size = 12, angle = 45, hjust = 1),
        plot.title = element_text(color = "black", size = 18, hjust = 0.5)) + 
  scale_color_manual(name = "Tissue type", labels = c("ND", "VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#D6A857", "#36DAFF", "#9C27B0", "#DF2E2E"))

dev.off()

# making boxplot with pattern 4 across the tissue types, with statistics 
newmetadatapat <- newmetadata %>% dplyr::select("VI_Classification", "Pattern_4")
data_for_p_val_manual <- p_table %>% filter(Pattern == "Pattern_4") %>% dplyr::select("group_1", "group_2", "p")
# data_for_p_val_manual$label <- paste0("p = ", signif(data_for_p_val_manual$p, digits = 3))
# data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "label")
colnames(data_for_p_val_manual) <- c("group1", "group2", "p")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)

# Convert p-values to asterisks
datman$p.signif <- sapply(datman$p, p_to_asterisk)

#datman$p <- signif(datman$p,digits=3)
yposition <- max(newmetadatapat[,"Pattern_4"]) * 1.05
stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
pdf("outputs/plots_tables_objects/Pattern_4_boxplot.pdf", width = 4, height = 6)
ggplot(newmetadata,aes(x = VI_Classification, y = Pattern_4, color = VI_Classification)) +
  geom_boxplot() + 
  geom_jitter() + 
  labs(color = "", x = "", y = "Pattern weight", title = "Pattern 4") +
  stat.test +
  theme(axis.text.x = element_text(color = "black",
                                   size = 12, angle = 45, hjust = 1),
        plot.title = element_text(color = "black", size = 18, hjust = 0.5)) + 
  scale_color_manual(name = "Tissue type", labels = c("ND", "VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#D6A857", "#36DAFF", "#9C27B0", "#DF2E2E"))
dev.off()
#####

# making boxplots of the pattern markers
jp4thresh <- readRDS('outputs/thresholded_patterns.rds')
annots <- read.csv("inputs/annotations.csv", header= F)
colnames(annots) <- c("Slide", "Number", "Classification")
annots <- as.data.frame(annots)
thedata <- read.csv("outputs/preprocessing/fully_batch_corrected_vsd.csv", row.names = 1)
metadata <- read.csv("outputs/preprocessing/metadata.csv", row.names = 1)
# removing stroma 063_003
annots <- annots[-87,]
metadata <- cbind(metadata, annots$Classification)
colnames(metadata)[16] <- "VI_Classification"
rownames(metadata) <- metadata$Name
newmetadata <- cbind(metadata,cogaps@sampleFactors)
newmetadata$VI_Classification[which(newmetadata$VI_Classification == "X")] <- "ND"
newmetadata$VI_Classification[which(newmetadata$VI_Classification == "A")] <- "VI-Destructive"
newmetadata$VI_Classification[which(newmetadata$VI_Classification == "B")] <- "VI-IN-Like"
newmetadata$VI_Classification[which(newmetadata$VI_Classification == "C")] <- "VI-Conventional"
newmetadata$VI_Classification[which(newmetadata$Type == "PDAC")] <- "PDAC"
newmetadata$VI_Classification[which(newmetadata$Type == "NORMAL")] <- "ND"
newmetadata$VI_Classification[which(newmetadata$Type == "PNI")] <- "PNI"

genesSig <- jp4thresh$PatternMarkers[[1]]
metadataDat <- data.frame(newmetadata,t(thedata[genesSig,]))
metadataDat <- metadataDat %>% filter(VI_Classification == "ND" | VI_Classification == "VI-Destructive" | VI_Classification == "VI-IN-Like" | VI_Classification == "PDAC")
metadataDat$VI_Classification <- factor(metadataDat$VI_Classification, levels = c("ND", "VI-IN-Like", "VI-Destructive", "PDAC"))

p_table_2 <- readRDS("outputs/plots_tables_objects/p_table_2.rds")

# FIGURE 6E, S8B-S8C ---- box plots of pattern marker genes between the VI subtypes ####
genesSig <- jp4thresh$PatternMarkers[[1]]
metadataDat <- data.frame(newmetadata,t(thedata[genesSig,]))
metadataDat <- metadataDat %>% filter(VI_Classification == "ND" | VI_Classification == "VI-Destructive" | VI_Classification == "VI-IN-Like" | VI_Classification == "PDAC")
metadataDat$VI_Classification <- factor(metadataDat$VI_Classification, levels = c("ND", "VI-IN-Like", "VI-Destructive", "PDAC"))
genesSig <- genesSig[which(genesSig %in% colnames(metadataDat))]
for (g in genesSig) {
  metadataDatgene <- metadataDat %>% dplyr::select("Type", g)
  data_for_p_val_manual <- p_table_2 %>% filter(Genes == make.names(g)) %>% dplyr::select("group_1", "group_2", "p", "LFC")
  #data_for_p_val_manual$label <- paste0("p.adj = ", signif(data_for_p_val_manual$p, digits = 2), ", ", "Log2FC = ", signif(data_for_p_val_manual$LFC, digits = 2))
  data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "p")
  colnames(data_for_p_val_manual) <- c("group1", "group2", "p.adj")
  data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
  datman <- tibble(data_for_p_val_manual)
  
  # Convert p-values to asterisks
  datman$p.signif <- sapply(datman$p.adj, p_to_asterisk)
  
  #datman$p.adj <- signif(datman$p.adj,digits=3)
  yposition <- max(metadataDatgene[,make.names(g)]) * 1.05
  stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
  pdf(paste0("outputs/geneplots/vi_subtypes/pattern_1_markers/", g,"_boxplot.pdf"), width = 4, height = 6)
  print(ggplot(metadataDat, aes_string(x='VI_Classification', y=make.names(g), color='VI_Classification')) + 
          geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
          labs(color = "", x = "", y = "Normalized expression", title = g) +
          theme(axis.text.x = element_text(color = "black",
                                           size = 12, angle = 45, hjust = 1),
                plot.title = element_text(hjust = 0.5, size = 18, face = "italic")) + 
          scale_color_manual(name = "Tissue type", labels = c("ND", "VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#D6A857", "#36DAFF", "#9C27B0", "#DF2E2E")) +
          stat.test)
  dev.off()
}

# making box plots for pattern 3 markers with stats
genesSig <- jp4thresh$PatternMarkers[[3]]
metadataDat <- data.frame(newmetadata,t(thedata[genesSig,]))
metadataDat <- metadataDat %>% filter(VI_Classification == "ND" | VI_Classification == "VI-Destructive" | VI_Classification == "VI-IN-Like" | VI_Classification == "PDAC")
metadataDat$VI_Classification <- factor(metadataDat$VI_Classification, levels = c("ND", "VI-IN-Like", "VI-Destructive", "PDAC"))
genesSig <- genesSig[which(genesSig %in% colnames(metadataDat))]
for (g in genesSig) {
  metadataDatgene <- metadataDat %>% dplyr::select("Type", g)
  data_for_p_val_manual <- p_table_2 %>% filter(Genes == make.names(g)) %>% dplyr::select("group_1", "group_2", "p", "LFC")
  #data_for_p_val_manual$label <- paste0("p.adj = ", signif(data_for_p_val_manual$p, digits = 2), ", ", "Log2FC = ", signif(data_for_p_val_manual$LFC, digits = 2))
  data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "p")
  colnames(data_for_p_val_manual) <- c("group1", "group2", "p.adj")
  data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
  datman <- tibble(data_for_p_val_manual)
  
  # Convert p-values to asterisks
  datman$p.signif <- sapply(datman$p.adj, p_to_asterisk)
  
  #datman$p.adj <- signif(datman$p.adj,digits=3)
  yposition <- max(metadataDatgene[,make.names(g)]) * 1.05
  stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
  pdf(paste0("outputs/geneplots/vi_subtypes/pattern_3_markers/", g,"_boxplot.pdf"), width = 4, height = 6)
  print(ggplot(metadataDat, aes_string(x='VI_Classification', y=make.names(g), color='VI_Classification')) + 
          geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
          labs(color = "", x = "", y = "Normalized expression", title = g) +
          theme(axis.text.x = element_text(color = "black",
                                           size = 12, angle = 45, hjust = 1),
                plot.title = element_text(hjust = 0.5, size = 18, face = "italic")) + 
          scale_color_manual(name = "Tissue type", labels = c("ND", "VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#D6A857", "#36DAFF", "#9C27B0", "#DF2E2E")) +
          stat.test)
  dev.off()
}

# making box plots for pattern 4 markers with stats
genesSig <- jp4thresh$PatternMarkers[[4]]
metadataDat <- data.frame(newmetadata,t(thedata[genesSig,]))
metadataDat <- metadataDat %>% filter(VI_Classification == "ND" | VI_Classification == "VI-Destructive" | VI_Classification == "VI-IN-Like" | VI_Classification == "PDAC")
metadataDat$VI_Classification <- factor(metadataDat$VI_Classification, levels = c("ND", "VI-IN-Like", "VI-Destructive", "PDAC"))
genesSig <- genesSig[which(genesSig %in% colnames(metadataDat))]
for (g in genesSig) {
  metadataDatgene <- metadataDat %>% dplyr::select("Type", g)
  data_for_p_val_manual <- p_table_2 %>% filter(Genes == make.names(g)) %>% dplyr::select("group_1", "group_2", "p")
  colnames(data_for_p_val_manual) <- c("group1", "group2", "p.adj")
  data_for_p_val_manual <- data_for_p_val_manual[c(1,3),]
  datman <- tibble(data_for_p_val_manual)
  
  # Convert p-values to asterisks
  datman$p.signif <- sapply(datman$p.adj, p_to_asterisk)
  
  datman$p.adj <- signif(datman$p.adj,digits=3)
  yposition <- max(metadataDatgene[,make.names(g)]) * 1.05
  stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
  pdf(paste0("outputs/geneplots/vi_subtypes/pattern_4_markers/", g,"_boxplot.pdf"), width = 4, height = 6)
  print(ggplot(metadataDat, aes_string(x='VI_Classification', y=make.names(g), color='VI_Classification')) + 
          geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
          labs(color = "", x = "", y = "Normalized expression", title = g) +
          theme(axis.text.x = element_text(color = "black",
                                           size = 12, angle = 45, hjust = 1),
                plot.title = element_text(hjust = 0.5, size = 18, face = "italic")) + 
          scale_color_manual(name = "Tissue type", labels = c("ND", "VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#D6A857", "#36DAFF", "#9C27B0", "#DF2E2E"))
        + stat.test)
  dev.off()
}

#####

# GSVA BETWEEN CLASSICAL AND BASAL BETWEEN THE DIFFERENT GROUPS (eg destructive, in-like, pdac)
# reading in batch corrected normalized data
thedata <- read.csv("outputs/preprocessing/fully_batch_corrected_vsd.csv", row.names = 1)
metadata <- read.csv("outputs/preprocessing/metadata.csv")

# reading in VI's annotation excel sheet
annots <- read.csv("inputs/annotations.csv", header= F)
colnames(annots) <- c("Slide", "Number", "Classification")
annots <- as.data.frame(annots)

# removing stroma 063_003
annots <- annots[-87,]

metadata <- cbind(metadata, annots$Classification)
colnames(metadata)[17] <- "VI_Classification"
rownames(metadata) <- metadata$Name

# reading in the cogaps object
cogaps <- readRDS("inputs/batch_merged_cogapsP4.rds")

# renaming the tissue subtypes
newmetadata <- cbind(metadata,cogaps@sampleFactors)
newmetadata$VI_Classification[which(newmetadata$VI_Classification == "A")] <- "VI-Destructive"
newmetadata$VI_Classification[which(newmetadata$VI_Classification == "B")] <- "VI-IN-Like"
newmetadata$VI_Classification[which(newmetadata$VI_Classification == "C")] <- "VI-conventional"
newmetadata$VI_Classification[which(newmetadata$Type == "PDAC")] <- "PDAC"
newmetadata$VI_Classification[which(newmetadata$Type == "NORMAL")] <- "Normal Duct"
newmetadata$VI_Classification[which(newmetadata$Type == "PNI")] <- "PNI"
newmetadata$VI_Classification

# getting only Destructive, IN-like, and PDAC
theABPindex <- which(newmetadata$VI_Classification == "VI-Destructive" | newmetadata$VI_Classification == "VI-IN-Like" | newmetadata$VI_Classification == "PDAC")

# filtering the data to only have destructive and IN-like samples
thedata <- thedata[,theABPindex]

# reading in the gene sets
collective <- c("TM4SF4", "ONECUT2", "RND1", "GATA6", "HNF1A")
mesenchymal <- c("ADAMTS9", "PTGS2", "TGFBI", "SERPINH1", "SPARC", "MFAP2", "MXRA8")

basal <- c('VGLL','UCA1','S100A2','LY6D','SPRR3','SPRR1B','LEMD1','KRT15',
           'CTSL2','DHRS9','AREG','CST6','SERPINB3','KRT6C','KRT6A','SERPINB4',
           'FAM83A','SCEL','FGFBP1','KRT7','KRT17','GPR87','TNS4','SLC2A1',
           'ANXA8L2')
classical <- c('BTNL8','FAM3D','ATAD4','AGR3','CTSE','LOC400573','LYZ','TFF2',
               'TFF1','ANXA10','LGALS4','PLA2G10','CEACAM6','VSIG2','TSPAN8',
               'ST6GALNAC1','AGR2','TFF3','CYP3A7','MYO1A','CLRN3','KRT20',
               'CDH17','SPINK4','REG4')
genesets <- list(basal, classical)
names(genesets) <- c("basal", "classical")

# performing GSVA for classical and basal-like
theoutput <- gsva(as.matrix(thedata), genesets, verbose=FALSE)

# merging the gsva results with the metadata
thetable <- t(theoutput)
themergednewmetadata <- merge(newmetadata, thetable, by = 0)

# reading the gsva results for specific tissue types
pdacbasal <-themergednewmetadata %>% filter(VI_Classification == "PDAC") %>% dplyr::select("basal") %>% as.vector() %>% unlist() %>% as.vector()
vdbasal <- themergednewmetadata%>% filter(VI_Classification == "VI-Destructive") %>% dplyr::select("basal") %>% as.vector() %>% unlist() %>% as.vector()
vinbasal <- themergednewmetadata %>% filter(VI_Classification == "VI-IN-Like") %>% dplyr::select("basal") %>% as.vector() %>% unlist() %>% as.vector()

pdacclassical <- themergednewmetadata %>% filter(VI_Classification == "PDAC") %>% dplyr::select("classical") %>% as.vector() %>% unlist() %>% as.vector()
vdclassical <- themergednewmetadata %>% filter(VI_Classification == "VI-Destructive") %>% dplyr::select("classical") %>% as.vector() %>% unlist() %>% as.vector()
vinclassical <- themergednewmetadata %>% filter(VI_Classification == "VI-IN-Like") %>% dplyr::select("classical") %>% as.vector() %>% unlist() %>% as.vector()


# doing the statistical tests between the tissue types
pdac_vdbasal <- wilcox.test(pdacbasal, vdbasal)
pdac_vinbasal <- wilcox.test(pdacbasal, vinbasal)
vin_vdbasal <- wilcox.test(vinbasal, vdbasal)

pdac_vdclassical <- wilcox.test(pdacclassical, vdclassical)
pdac_vinclassical <- wilcox.test(pdacclassical, vinclassical)
vin_vdclassical <- wilcox.test(vinclassical, vdclassical)

# making data frames of the p values
pdac_vdbasalp <- data.frame("Pattern" = "basal", "group_1" = "PDAC", "group_2" = "VI-Destructive", "p" = pdac_vdbasal$p.value)
pdac_vinbasalp <- data.frame("Pattern" = "basal", "group_1" = "PDAC", "group_2" = "VI-IN-Like", "p" = pdac_vinbasal$p.value)
vin_vdbasalp <- data.frame("Pattern" = "basal", "group_1" = "VI-IN-Like", "group_2" = "VI-Destructive", "p" = vin_vdbasal$p.value)

pdac_vdclassicalp <- data.frame("Pattern" = "classical", "group_1" = "PDAC", "group_2" = "VI-Destructive", "p" = pdac_vdclassical$p.value)
pdac_vinclassicalp <- data.frame("Pattern" = "classical", "group_1" = "PDAC", "group_2" = "VI-IN-Like", "p" = pdac_vinclassical$p.value)
vin_vdclassicalp <- data.frame("Pattern" = "classical", "group_1" = "VI-IN-Like", "group_2" = "VI-Destructive", "p" = vin_vdclassical$p.value)

# making a table of p values
p_table <- rbind(pdac_vdbasalp, pdac_vinbasalp, vin_vdbasalp, pdac_vdclassicalp, pdac_vinclassicalp, vin_vdclassicalp)

# FIGURE S6D ---- making boxplots of GSVA results for classical and basal-like ####
# making a boxplot of GSVA results for basal-like
newmetadatapat <- themergednewmetadata %>% dplyr::select("VI_Classification", "basal")
data_for_p_val_manual <- p_table %>% filter(Pattern == "basal") %>% dplyr::select("group_1", "group_2", "p")
#data_for_p_val_manual$label <- paste0("p = ", signif(data_for_p_val_manual$p, digits = 3))
data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "p")
colnames(data_for_p_val_manual) <- c("group1", "group2", "p")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)

# Convert p-values to asterisks
datman$p.signif <- sapply(datman$p, p_to_asterisk)

#datman$p <- signif(datman$p,digits=3)
yposition <- max(newmetadatapat[,"basal"]) * 1.05
stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
pdf("outputs/plots_tables_objects/basal_VI_gsva.pdf", width = 4, height = 6)
ggplot(themergednewmetadata, aes(x = VI_Classification, y = basal, color = VI_Classification)) + 
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(color = "", x = "", y = "GSVA enrichment score") +
  stat.test +
  theme(axis.text.x = element_text(color = "black",
                                   size = 12, angle = 45, hjust = 1)) + 
  scale_color_manual(name = "Tissue type", labels = c("PDAC", "VI-Destructive", "VI-IN-Like"), values = c("#DF2E2E", "#9C27B0", "#36DAFF"))
dev.off()

# making a boxplot of GSVA results for classical 
newmetadatapat <- themergednewmetadata %>% dplyr::select("VI_Classification", "classical")
data_for_p_val_manual <- p_table %>% filter(Pattern == "classical") %>% dplyr::select("group_1", "group_2", "p")
#data_for_p_val_manual$label <- paste0("p = ", signif(data_for_p_val_manual$p, digits = 3))
data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "p")
colnames(data_for_p_val_manual) <- c("group1", "group2", "p")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)

# Convert p-values to asterisks
datman$p.signif <- sapply(datman$p, p_to_asterisk)

#datman$p <- signif(datman$p,digits=3)
yposition <- max(newmetadatapat[,"classical"]) * 1.05
stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
pdf("outputs/plots_tables_objects/classical_VI_gsva.pdf", width = 4, height = 6)
ggplot(themergednewmetadata, aes(x = VI_Classification, y = classical, color = VI_Classification)) + 
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(color = "", x = "", y = "GSVA enrichment score") +
  stat.test +
  theme(axis.text.x = element_text(color = "black",
                                   size = 12, angle = 45, hjust = 1)) + 
  scale_color_manual(name = "Tissue type", labels = c("PDAC", "VI-Destructive", "VI-IN-Like"), values = c("#DF2E2E", "#9C27B0", "#36DAFF"))
dev.off()

#####

# now preparing to do the GSVA with mesenchymal and collective subtypes
genesets <- list(mesenchymal, collective)
names(genesets) <- c("mesenchymal", "collective")

# doing GSVA with mesenchymal and collective, and merging it with the metadata 
theoutput <- gsva(as.matrix(thedata), genesets, verbose=FALSE)
thetable <- t(theoutput)
themergednewmetadata <- merge(newmetadata, thetable, by = 0)

# isolating the mesenchymal GSVA scores for each subtype
pdacmesenchymal <- themergednewmetadata %>% filter(VI_Classification == "PDAC") %>% dplyr::select("mesenchymal") %>% as.vector() %>% unlist() %>% as.vector()
vdmesenchymal <- themergednewmetadata %>% filter(VI_Classification == "VI-Destructive") %>% dplyr::select("mesenchymal") %>% as.vector() %>% unlist() %>% as.vector()
vinmesenchymal <- themergednewmetadata %>% filter(VI_Classification == "VI-IN-Like") %>% dplyr::select("mesenchymal") %>% as.vector() %>% unlist() %>% as.vector()

# isolating the collective GSVA scores for each subtype
pdaccollective <- themergednewmetadata %>% filter(VI_Classification == "PDAC") %>% dplyr::select("collective") %>% as.vector() %>% unlist() %>% as.vector()
vdcollective <- themergednewmetadata %>% filter(VI_Classification == "VI-Destructive") %>% dplyr::select("collective") %>% as.vector() %>% unlist() %>% as.vector()
vincollective <- themergednewmetadata %>% filter(VI_Classification == "VI-IN-Like") %>% dplyr::select("collective") %>% as.vector() %>% unlist() %>% as.vector()

# doing stats for mesenchymal
pdac_vdmesenchymal <- wilcox.test(pdacmesenchymal, vdmesenchymal)
pdac_vinmesenchymal <- wilcox.test(pdacmesenchymal, vinmesenchymal)
vin_vdmesenchymal <- wilcox.test(vinmesenchymal, vdmesenchymal)

# doing stats for collective
pdac_vdcollective <- wilcox.test(pdaccollective, vdcollective)
pdac_vincollective <- wilcox.test(pdaccollective, vincollective)
vin_vdcollective <- wilcox.test(vincollective, vdcollective)

# making p value data frames for mesenchymal
pdac_vdmesenchymalp <- data.frame("Pattern" = "mesenchymal", "group_1" = "PDAC", "group_2" = "VI-Destructive", "p" = pdac_vdmesenchymal$p.value)
pdac_vinmesenchymalp <- data.frame("Pattern" = "mesenchymal", "group_1" = "PDAC", "group_2" = "VI-IN-Like", "p" = pdac_vinmesenchymal$p.value)
vin_vdmesenchymalp <- data.frame("Pattern" = "mesenchymal", "group_1" = "VI-IN-Like", "group_2" = "VI-Destructive", "p" = vin_vdmesenchymal$p.value)

# making p value data frames for collective
pdac_vdcollectivep <- data.frame("Pattern" = "collective", "group_1" = "PDAC", "group_2" = "VI-Destructive", "p" = pdac_vdcollective$p.value)
pdac_vincollectivep <- data.frame("Pattern" = "collective", "group_1" = "PDAC", "group_2" = "VI-IN-Like", "p" = pdac_vincollective$p.value)
vin_vdcollectivep <- data.frame("Pattern" = "collective", "group_1" = "VI-IN-Like", "group_2" = "VI-Destructive", "p" = vin_vdcollective$p.value)

# combining all the p value tables
p_table <- rbind(pdac_vdmesenchymalp, pdac_vinmesenchymalp, vin_vdmesenchymalp, pdac_vdcollectivep, pdac_vincollectivep, vin_vdcollectivep)

# FIGURE S6E ---- # making boxplot for collective and mesenchymal signature across all VI subtypes #####
# collective signature boxplot
newmetadatapat <- themergednewmetadata %>% dplyr::select("VI_Classification", "collective")
data_for_p_val_manual <- p_table %>% filter(Pattern == "collective") %>% dplyr::select("group_1", "group_2", "p")
#data_for_p_val_manual$label <- paste0("p = ", signif(data_for_p_val_manual$p, digits = 3))
data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "p")
colnames(data_for_p_val_manual) <- c("group1", "group2", "p")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)

# Convert p-values to asterisks
datman$p.signif <- sapply(datman$p, p_to_asterisk)

#datman$p <- signif(datman$p,digits=3)
yposition <- max(newmetadatapat[,"collective"]) * 1.05
stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
pdf("outputs/plots_tables_objects/collective_VI_gsva.pdf", width = 4, height = 6)
ggplot(themergednewmetadata, aes(x = VI_Classification, y = collective, color = VI_Classification)) + 
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(color = "", x = "", y = "GSVA enrichment score") +
  stat.test +
  theme(axis.text.x = element_text(color = "black",
                                   size = 12, angle = 45, hjust = 1)) +
  scale_color_manual(name = "Tissue type", labels = c("PDAC", "VI-Destructive", "VI-IN-Like"), values = c("#DF2E2E", "#9C27B0", "#36DAFF"))
dev.off()

# mesenchymal signature boxplot
newmetadatapat <- themergednewmetadata %>% dplyr::select("VI_Classification", "mesenchymal")
data_for_p_val_manual <- p_table %>% filter(Pattern == "mesenchymal") %>% dplyr::select("group_1", "group_2", "p")
#data_for_p_val_manual$label <- paste0("p = ", signif(data_for_p_val_manual$p, digits = 3))
data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "p")
colnames(data_for_p_val_manual) <- c("group1", "group2", "p")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)

# Convert p-values to asterisks
datman$p.signif <- sapply(datman$p, p_to_asterisk)

#datman$p <- signif(datman$p,digits=3)
yposition <- max(newmetadatapat[,"mesenchymal"]) * 1.05
stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
pdf("outputs/plots_tables_objects/mesenchymal_VI_gsva.pdf", width = 4, height = 6)
ggplot(themergednewmetadata, aes(x = VI_Classification, y = mesenchymal, color = VI_Classification)) + 
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
  labs(color = "", x = "", y = "GSVA enrichment score") +
  stat.test +
  theme(axis.text.x = element_text(color = "black",
                                   size = 12, angle = 45, hjust = 1)) +
  scale_color_manual(name = "Tissue type", labels = c("PDAC", "VI-Destructive", "VI-IN-Like"), values = c("#DF2E2E", "#9C27B0", "#36DAFF"))
dev.off()

#####

# Extended snRNA Seq phenotypes
# Reading in the snRNAseq gene sets
cell_types <- read.csv("inputs/snRNA_seq_cell_type_sets.csv")
cell_states <- read.csv("inputs/snRNA_seq_cell_state_sets.csv")

# cell_types[,0:7] because some reason there are two extra columns being imported
cell_types_list <- lapply(cell_types[,0:7], as.vector)
cell_states_list <- lapply(cell_states, as.vector)

# reading in the gene sets
genesets <- cell_types_list

# performing GSVA for cell types
theoutput <- gsva(as.matrix(thedata), genesets, verbose=FALSE)

# merging the gsva results with the metadata
thetable <- t(theoutput)
themergednewmetadata <- merge(newmetadata, thetable, by = 0)

# FIGURE S7A ---- making boxplots of GSVA results for the Hwang et al 2022 cell lineages ####
for(pathway in names(genesets)){
  pdf(paste0("outputs/plots_tables_objects/snRNAseq/",pathway,".pdf"), width = 4, height = 6)
  print(ggplot(themergednewmetadata, aes_string(x = "VI_Classification", y = paste0(pathway), color = "VI_Classification")) +
          geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
          labs(color = "", x = "", y = "GSVA enrichment score", title = pathway) +
          # stat.test +
          theme(axis.text.x = element_text(color = "black",
                                           size = 12, angle = 45, hjust = 1)) +
          scale_color_manual(name = "Tissue type", labels = c("PDAC", "VI-Destructive", "VI-IN-Like"), values = c("#DF2E2E", "#9C27B0", "#36DAFF")) +
          stat_compare_means(comparisons = list(
            c("PDAC", "VI-Destructive"),
            c("PDAC", "VI-IN-Like"),
            c("VI-Destructive", "VI-IN-Like")), method = "wilcox.test", label = "p.signif", size = 4))
  dev.off()
}

# reading in the gene sets
genesets <- cell_states_list

# performing GSVA for cell types
theoutput <- gsva(as.matrix(thedata), genesets, verbose=FALSE)

# merging the gsva results with the metadata
thetable <- t(theoutput)
themergednewmetadata <- merge(newmetadata, thetable, by = 0)

# FIGURE S7B ---- making boxplots of GSVA results for the Hwang et al cell states ####
for(pathway in names(genesets)){
  pdf(paste0("outputs/plots_tables_objects/snRNAseq/",pathway,".pdf"), width = 4, height = 6)
  print(ggplot(themergednewmetadata, aes_string(x = "VI_Classification", y = paste0(pathway), color = "VI_Classification")) +
          geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
          labs(color = "", x = "", y = "GSVA enrichment score", title = pathway) +
          # stat.test +
          theme(axis.text.x = element_text(color = "black",
                                           size = 12, angle = 45, hjust = 1)) +
          scale_color_manual(name = "Tissue type", labels = c("PDAC", "VI-Destructive", "VI-IN-Like"), values = c("#DF2E2E", "#9C27B0", "#36DAFF")) +
          stat_compare_means(comparisons = list(
            c("PDAC", "VI-Destructive"),
            c("PDAC", "VI-IN-Like"),
            c("VI-Destructive", "VI-IN-Like")), method = "wilcox.test", label = "p.signif", size = 4))
  dev.off()
}

#####

