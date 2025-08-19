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
library(CoGAPS)
library(fgsea)
library(msigdbr)
library(GSVA)
library("reshape2")
library(scales)  
sessionInfo()

source("scripts/00_Custom_Functions.R")
##### Importing data, outputs from script 01
batch_adjusted <- readRDS("outputs/preprocessing/batch_adjusted.rds")
metadata <- readRDS("outputs/preprocessing/metadata.rds")
pheno = new("AnnotatedDataFrame", data=metadata)
eset = ExpressionSet(assayData = batch_adjusted, phenoData = pheno)

# Differential expression analysis between conditions

# Making matrix design
design = model.matrix(~0 + pheno$Type)
colnames(design) <- levels(pheno$Type)

# Making contrasts for differential expression
contr.matrix <- makeContrasts(
  VIvsND = VI - ND,
  VIvsPNI = VI - PNI, 
  VIvsPDAC = VI - PDAC,
  PNIvsND = PNI - ND,
  PNIvsPDAC= PNI - PDAC,
  PDACvsND = PDAC - ND,
  levels = colnames(design))

# doing limma differential expression on the batch adjusted matrix
eset = ExpressionSet(assayData = batch_adjusted, phenoData = pheno)
fit <- lmFit(eset, design)
vfit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(vfit)

# topTreat 
vi.vs.normal <- topTreat(efit, coef=1, n=Inf)
vi.vs.pni <- topTreat(efit, coef = 2, n =Inf)
vi.vs.pdac <- topTreat(efit, coef=3, n=Inf)
saveRDS(vi.vs.pdac, "outputs/plots_tables_objects/vipdaclimma.rds")
pni.vs.normal <- topTreat(efit, coef = 4, n =Inf)
pni.vs.pdac <- topTreat(efit, coef = 5, n =Inf)
pdac.vs.normal <- topTreat(efit, coef = 6, n =Inf)

# FIGURE 3A ---- Limma differential expression analysis between VI and PDAC samples ####
volcanoPlot(vi.vs.pdac, title_1 = "VI", title_2 = "PDAC", upcolor = "#2196F3", downcolor = "#DF2E2E", batch = "Differential_Expression", directory = "outputs/plots_tables_objects/")
#####

# make other volcano plots
volcanoPlot(vi.vs.normal, title_1 = "VI", title_2 = "ND", batch = "Differential_Expression", directory = "outputs/random/")
volcanoPlot(vi.vs.pni, title_1 = "VI", title_2 = "PNI", batch = "Differential_Expression", directory = "outputs/random/")
volcanoPlot(pni.vs.normal, title_1 = "PNI", title_2 = "ND", batch = "Differential_Expression", directory = "outputs/random/")
volcanoPlot(pni.vs.pdac, title_1 = "PNI", title_2 = "PDAC", batch = "Differential_Expression", directory = "outputs/random/")
volcanoPlot(pdac.vs.normal, title_1 = "PDAC", title_2 = "ND", batch = "Differential_Expression", directory = "outputs/random/")

# Getting adjusted p values for each comparison -- this will be used later in the creation of boxplots

vi_nd <- data.frame("Genes" = rownames(vi.vs.normal), "group_1" = "VI", "group_2" = "ND", "p" = vi.vs.normal$adj.P.Val)
vi_pni <- data.frame("Genes" = rownames(vi.vs.pni), "group_1" = "VI", "group_2" = "PNI", "p" = vi.vs.pni$adj.P.Val)
vi_pdac <- data.frame("Genes" = rownames(vi.vs.pdac), "group_1" = "VI", "group_2" = "PDAC", "p" = vi.vs.pdac$adj.P.Val)
pni_nd <- data.frame("Genes" = rownames(pni.vs.normal), "group_1" = "PNI", "group_2" = "ND", "p" = pni.vs.normal$adj.P.Val)
pni_pdac <- data.frame("Genes" = rownames(pni.vs.pdac), "group_1" = "PNI", "group_2" = "PDAC", "p" = pni.vs.pdac$adj.P.Val)
pdac_nd <- data.frame("Genes" = rownames(pdac.vs.normal), "group_1" = "PDAC", "group_2" = "ND", "p" = pdac.vs.normal$adj.P.Val)

# combining p values
p_table <- rbind(vi_nd, vi_pni, vi_pdac, pni_nd, pni_pdac, pdac_nd)
saveRDS(p_table, "outputs/plots_tables_objects/p_table.rds")

# Making volcano plots comparing all of the cancers to the normals
# Making matrix design
design = model.matrix(~0 + pheno$Cancer)
colnames(design) <- levels(pheno$Cancer)
# Making contrasts for differential expression
contr.matrix <- makeContrasts(
  PDACvsNormal = COMBINED_PDAC - ND,
  levels = colnames(design))
# doing limma differential expression on the batch adjusted matrix
eset = ExpressionSet(assayData = batch_adjusted, phenoData = pheno)
fit <- lmFit(eset, design)
vfit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(vfit)

# combined_pdac vs normal
cancer.vs.normal <- topTreat(efit, coef=1, n=Inf)

# FIGURE 2C ---- Differential expression analysis between Combined PDAC samples (PNI + PDAC + VI) and normal duct (ND). #####
volcanoPlot(cancer.vs.normal, title_1 = "Combined PDAC", title_2 = "ND", batch = "Differential_Expression", directory = "outputs/plots_tables_objects/")
#####

#filtering out significant genes (used 0.58 as lfc cutoff for volcano plots)
significant_degs <- vi.vs.pdac %>% filter(adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))

# making batch adjusted table of significant degs
significant_degs_table <- batch_adjusted[rownames(batch_adjusted) %in% rownames(significant_degs),]

# prepping data for heatmap # note that am just doing ALL genes instead of significant degs
datScale <- t(apply(batch_adjusted,1,scale))
colnames(datScale) <- colnames(batch_adjusted)
column_order = order(metadata$Cancer)

# These are the genes that will appear highlighted on the heatmap
thegenes <- c("FXYD2", "CFTR", "SLC3A1", "SLC4A4", "ANPEP", "CLU", "GPRC5A", "S100P", "IFI27", "FXYD3", "TSPAN1", "TMPRSS4")
theindex <- which(rownames(datScale) %in% thegenes)
thegenes <- rownames(datScale)[theindex]

anno = anno_mark(at = theindex, labels = thegenes, which = "row",
                 labels_gp = gpar(fontface = "italic"))
cancervsnormaldegs <- read.csv("outputs/plots_tables_objects/Differential_Expression_Combined PDACvsND_all_genes.csv")
rownames(cancervsnormaldegs) <- cancervsnormaldegs$rn
groupindex <- match(rownames(datScale), rownames(cancervsnormaldegs))
groupannonames <- cancervsnormaldegs$Group[groupindex]
row_ha = rowAnnotation(DEGs = groupannonames, col = list(DEGs = c("Uppers" = "#619CFF", "Mids" = "#F0FFF0", "Lowers" = "#F8766D")),
                       annotation_legend_param = list(DEGs = 
                                                        list(title = "",
                                                             at = c("Lowers", "Mids", "Uppers"),
                                                             labels = c("Increased in ND", "Unchanged", "Increased in Combined PDAC"))))

# FIGURE 2D ---- making heatmap of the differentially expressed genes between all cancer and normal duct ####
pdf(file = "outputs/plots_tables_objects/CombinedCancer_vs_normal_DEGS_heatmap_full.pdf")
Heatmap(datScale, show_row_names = F, show_column_names = F,
        name = "Matrix",
        clustering_distance_rows = 'pearson',
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = HeatmapAnnotation(Type =
                                             metadata[,c('Type')],
                                           Patient.Alias = metadata[,c('Patient.Alias')], col = list(Type= c("ND" = "#D6A857", "PDAC" = "#DF2E2E", "PNI" = "#E475FF", "VI" = "#2196F3"),
                                                                                               Patient.Alias = c("1" = "#00BCD4", "2" = "#000000", "3" = "#CBCBCB", "4" = "#9F5500", "5" = "#7AFF88", "6" = "#3F51B5", "7" = "#FF6600", "8" = "#FFBAE2")),
                                           annotation_legend_param = list(Patient.Alias = list(title = "Patient alias"), Type = list(title = "Tissue type"))), 
        right_annotation = row_ha, 
        heatmap_legend_param = list()) + rowAnnotation(mark = anno)
dev.off()
#####


# Plot heatmap of VI vs PDAC genes
#filtering out significant genes (used 0.58 as lfc cutoff for volcano plots)
significant_degs <- vi.vs.pdac %>% filter(adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58))

# making batch adjusted table of significant degs
significant_degs_table <- batch_adjusted[rownames(batch_adjusted) %in% rownames(significant_degs),]

# prepping data for heatmap
datScale <- t(apply(significant_degs_table,1,scale))
colnames(datScale) <- colnames(batch_adjusted)
column_order = order(metadata$Type)

# these are the list of genes that get automatically labeled on the volcano plot of VI vs PDAC.. to be row annotation for the heatmap
thegenes <- c("MMP14", "PLP2", "LUM", "PLAU", "MMP11", "TGM2", "FN1", "CYP2C9", "HSD17B2", "MUC13", "CEACAM5", "CFTR", "REG1A", "PSCA")

# reading in the DEG table of VI vs PDAC
vivspdacdegs <- read.csv("outputs/plots_tables_objects/Differential_Expression_VIvsPDAC_all_genes.csv")

# making row annotation for the heatmap
theindex <- which(rownames(datScale) %in% thegenes)
thegenes <- rownames(datScale)[theindex]
anno = anno_mark(at = theindex, labels = thegenes, which = "row",
                 labels_gp = gpar(fontface = "italic"))
rownames(vivspdacdegs) <- vivspdacdegs$rn
groupindex <- match(rownames(datScale), rownames(vivspdacdegs))
groupannonames <- vivspdacdegs$Group[groupindex]
row_ha = rowAnnotation(DEGs = groupannonames, col = list(DEGs = c("Uppers" = "#619CFF", "Lowers" = "#F8766D")),
                       annotation_legend_param = list(DEGs = 
                                                        list(title = "",
                                                             at = c("Lowers", "Uppers"),
                                                             labels = c("Increased in PDAC", "Increased in VI"))))

# FIGURE 3B ---- Heatmap of samples grouped by tissue type, and significant VI-PDAC DEGs clustered by Pearson correlation #### 
# making heatmap of VI vs PDAC with row annotation of the genes labeled on the volcano plot
pdf(file = "outputs/plots_tables_objects/VIPDAC_DEGS_heatmap_full.pdf")
Heatmap(datScale, show_row_names = F, show_column_names = F,
        name = "Matrix",
        clustering_distance_rows = 'pearson',
        column_order = column_order,
        row_names_gp = grid::gpar(fontsize = 8),
        top_annotation = HeatmapAnnotation(Type =
                                             metadata[,c('Type')],
                                           Patient.Alias = metadata[,c('Patient.Alias')], col = list(Type= c("ND" = "#D6A857", "PDAC" = "#DF2E2E", "PNI" = "#E475FF", "VI" = "#2196F3"),
                                                                                               Patient.Alias = c("1" = "#00BCD4", "2" = "#000000", "3" = "#CBCBCB", "4" = "#9F5500", "5" = "#7AFF88", "6" = "#3F51B5", "7" = "#FF6600", "8" = "#FFBAE2")),
                                           annotation_legend_param = list(Patient.Alias = list(title = "Patient alias"), Type = list(title = "Tissue type"))), 
        right_annotation = row_ha) + rowAnnotation(mark = anno)
dev.off()

#####

# RUNNING COGAPS

# geoData <- as.matrix(outputs/preprocessing/fully_batch_corrected_vsd.csv", row.names = 1))
# geoGapsResultsP4 <- CoGAPS(data = geoData, nPatterns = 4, nIterations = 50000)
# saveRDS(geoGapsResultsP4,"outputs/plots_tables_objects/batch_merged_cogapsp4.csv")

## MAKING COGAPS BOX PLOTS

# reading in CoGAPS results
geoGapsResultsP4 <- readRDS(file = "inputs/batch_merged_cogapsP4.rds")
geoGapsResultsP4@sampleFactors
jp4thresh <- patternMarkers(geoGapsResultsP4, threshold="cut")
patternmarkerz <- t(plyr::ldply(jp4thresh$PatternMarkers, rbind))
saveRDS(jp4thresh, "outputs/thresholded_patterns.rds")

colnames(patternmarkerz) <- c("Pattern 1", "Pattern 2", "Pattern 3", "Pattern 4")
write.csv(geoGapsResultsP4@sampleFactors, "outputs/sampleFactors.csv")
write.csv(geoGapsResultsP4@featureLoadings, "outputs/featureLoadings.csv")
write.csv(patternmarkerz, "outputs/pattern_markers.csv")

tissue_type <- c("VI", "VI", "ND",
                 "VI", "VI", "ND",
                 "PDAC", "PDAC", "PDAC",
                 "PDAC", "PDAC", "PDAC",
                 "VI", "VI", "VI", 
                 "ND", "ND", "ND",
                 "ND", "PDAC", "PDAC",
                 "PDAC", "PDAC", "PDAC",
                 "VI", "VI", "VI",
                 "VI", "VI", "PDAC",
                 "ND", "ND", "PDAC",
                 "PDAC", "PDAC", "PNI",
                 "VI", "PNI", "PNI",
                 "VI", "VI", "VI",
                 "ND", "ND", "PDAC",
                 "PDAC", "PDAC", "PDAC",
                 
                 "VI", "VI", "VI",
                 "VI", "VI", "VI",
                 "PNI", "VI", "PDAC",
                 "PDAC", "PDAC", "PDAC",
                 "VI", "ND", "ND",
                 "ND", "VI", "VI",
                 "PDAC", "PDAC", "PDAC",
                 "PDAC", "PDAC", "PDAC",
                 "VI", "VI", "PDAC",
                 "PDAC", "PDAC", "PDAC",
                 "PDAC", "PDAC", "PDAC",
                 "VI", "VI", "VI",
                 "VI", "PNI", "PNI", 
                 "PNI", "PNI", "VI",
                 "VI", "PDAC", "PDAC",
                 "VI", "ND")

# Slide for merged   
slide <- c("1", "1", "1",
           "1", "1", "1",
           "1", "1", "1",
           "1", "1", "1",
           "2", "2", "2",
           "2", "2", "2",
           "2", "2", "2",
           "2", "2", "2",
           "3", "3", "3",
           "3", "3", "3",
           "3", "3", "3",
           "3", "3", "3",
           "4", "4", "4",
           "4", "4", "4",
           "4", "4", "4",
           "4", "4", "4",
           "5", "5", "5",
           "5", "5", "5",
           "5", "5", "5",
           "5", "5", "5",
           "6", "6", "6",
           "6", "6", "6",
           "6", "6", "6",
           "6", "6", "6",
           "7", "7", "7",
           "7", "7", "7",
           "7", "7", "7",
           "7", "7", "7",
           "8", "8", "8",
           "8", "8", "8",
           "8", "8", "8",
           "8", "8")


tempDf <- as.data.frame(geoGapsResultsP4@sampleFactors)
# make table of tissue type and slide
type_and_slide <- cbind(tissue_type, slide)
# combine metadata with pattern matrix
colDataDf <- cbind(tempDf, type_and_slide)
#reordering the columns from the boxplot
colDataDf$tissue_type <- factor(colDataDf$tissue_type, levels=c('ND','PDAC', 'PNI','VI'))


# this block of code is isolating the patterns by tissue type, for later use in pairwise statistical comparisons
vi1 <- colDataDf %>% filter(tissue_type == "VI") %>% dplyr::select("Pattern_1") %>% as.vector() %>% unlist() %>% as.vector()
pdac1 <- colDataDf %>% filter(tissue_type == "PDAC") %>% dplyr::select("Pattern_1") %>% as.vector() %>% unlist() %>% as.vector()
pni1 <- colDataDf %>% filter(tissue_type == "PNI") %>% dplyr::select("Pattern_1") %>% as.vector() %>% unlist() %>% as.vector()
nd1 <- colDataDf %>% filter(tissue_type == "ND") %>% dplyr::select("Pattern_1") %>% as.vector() %>% unlist() %>% as.vector()
vi2 <- colDataDf %>% filter(tissue_type == "VI") %>% dplyr::select("Pattern_2") %>% as.vector() %>% unlist() %>% as.vector()
pdac2 <- colDataDf %>% filter(tissue_type == "PDAC") %>% dplyr::select("Pattern_2") %>% as.vector() %>% unlist() %>% as.vector()
pni2 <- colDataDf %>% filter(tissue_type == "PNI") %>% dplyr::select("Pattern_2") %>% as.vector() %>% unlist() %>% as.vector()
nd2 <- colDataDf %>% filter(tissue_type == "ND") %>% dplyr::select("Pattern_2") %>% as.vector() %>% unlist() %>% as.vector()
vi3 <- colDataDf %>% filter(tissue_type == "VI") %>% dplyr::select("Pattern_3") %>% as.vector() %>% unlist() %>% as.vector()
pdac3 <- colDataDf %>% filter(tissue_type == "PDAC") %>% dplyr::select("Pattern_3") %>% as.vector() %>% unlist() %>% as.vector()
pni3 <- colDataDf %>% filter(tissue_type == "PNI") %>% dplyr::select("Pattern_3") %>% as.vector() %>% unlist() %>% as.vector()
nd3 <- colDataDf %>% filter(tissue_type == "ND") %>% dplyr::select("Pattern_3") %>% as.vector() %>% unlist() %>% as.vector()
vi4 <- colDataDf %>% filter(tissue_type == "VI") %>% dplyr::select("Pattern_4") %>% as.vector() %>% unlist() %>% as.vector()
pdac4 <- colDataDf %>% filter(tissue_type == "PDAC") %>% dplyr::select("Pattern_4") %>% as.vector() %>% unlist() %>% as.vector()
pni4 <- colDataDf %>% filter(tissue_type == "PNI") %>% dplyr::select("Pattern_4") %>% as.vector() %>% unlist() %>% as.vector()
nd4 <- colDataDf %>% filter(tissue_type == "ND") %>% dplyr::select("Pattern_4") %>% as.vector() %>% unlist() %>% as.vector()


# this block of code is performing wilcoxon ranked sum tests between each tissue type of interest for each cogaps pattern
p_vipdac1 <- wilcox.test(vi1, pdac1)
p_vipdac1 <- p_vipdac1$p.value
p_vipdac1 <- data.frame("Pattern" = "1", "group_1" = "VI", "group_2" = "PDAC", "p" = p_vipdac1)
p_vind1 <- wilcox.test(vi1, nd1)
p_vind1 <- p_vind1$p.value
p_vind1 <- data.frame("Pattern" = "1", "group_1" = "VI", "group_2" = "ND", "p" = p_vind1)
p_pdacnd1 <- wilcox.test(pdac1, nd1)
p_pdacnd1 <- p_pdacnd1$p.value
p_pdacnd1 <- data.frame("Pattern" = "1", "group_1" = "PDAC", "group_2" = "ND", "p" = p_pdacnd1)
p_vipdac2 <- wilcox.test(vi2, pdac2)
p_vipdac2 <- p_vipdac2$p.value
p_vipdac2 <- data.frame("Pattern" = "2", "group_1" = "VI", "group_2" = "PDAC", "p" = p_vipdac2)
p_vind2 <- wilcox.test(vi2, nd2)
p_vind2 <- p_vind2$p.value
p_vind2 <- data.frame("Pattern" = "2", "group_1" = "VI", "group_2" = "ND", "p" = p_vind2)
p_pdacnd2 <- wilcox.test(pdac2, nd2)
p_pdacnd2 <- p_pdacnd2$p.value
p_pdacnd2 <- data.frame("Pattern" = "2", "group_1" = "PDAC", "group_2" = "ND", "p" = p_pdacnd2)
p_vipdac3 <- wilcox.test(vi3, pdac3)
p_vipdac3 <- p_vipdac3$p.value
p_vipdac3 <- data.frame("Pattern" = "3", "group_1" = "VI", "group_2" = "PDAC", "p" = p_vipdac3)
p_vind3 <- wilcox.test(vi3, nd3)
p_vind3 <- p_vind3$p.value
p_vind3 <- data.frame("Pattern" = "3", "group_1" = "VI", "group_2" = "ND", "p" = p_vind3)
p_pdacnd3 <- wilcox.test(pdac3, nd3)
p_pdacnd3 <- p_pdacnd3$p.value
p_pdacnd3 <- data.frame("Pattern" = "3", "group_1" = "PDAC", "group_2" = "ND", "p" = p_pdacnd3)
p_vipdac4 <- wilcox.test(vi4, pdac4)
p_vipdac4 <- p_vipdac4$p.value
p_vipdac4 <- data.frame("Pattern" = "4", "group_1" = "VI", "group_2" = "PDAC", "p" = p_vipdac4)
p_vind4 <- wilcox.test(vi4, nd4)
p_vind4 <- p_vind4$p.value
p_vind4 <- data.frame("Pattern" = "4", "group_1" = "VI", "group_2" = "ND", "p" = p_vind4)
p_pdacnd4 <- wilcox.test(pdac4, nd4)
p_pdacnd4 <- p_pdacnd4$p.value
p_pdacnd4 <- data.frame("Pattern" = "4", "group_1" = "PDAC", "group_2" = "ND", "p" = p_pdacnd4)

# combining all the p values into one table, for later use in boxplot annotation
p_table <- rbind(p_vind1, p_vind2, p_vind3, p_vind4, p_vipdac1, p_vipdac2, p_vipdac3, p_vipdac4, p_pdacnd1, p_pdacnd2, p_pdacnd3, p_pdacnd4)

# FIGURE 3C ---- Boxplots of each of the four gene expression patterns discovered by CoGAPS from the normalized expression matrix. Each pattern is plotted by weight per tissue type.####
# creating boxplots of each cogaps pattern by tissue type
patterns <- colnames(tempDf)
for(x in 1:length(patterns)){
  data_for_p_val_manual <- p_table %>% filter(Pattern == x) %>% dplyr::select("group_1", "group_2", "p")
  colnames(data_for_p_val_manual) <- c("group1", "group2", "p")
  datman <- tibble(data_for_p_val_manual)
  
  # Convert p-values to asterisks
  datman$p.signif <- sapply(datman$p, p_to_asterisk)
  
  datman$p <- signif(datman$p,digits=3)
  yposition <- max(colDataDf[,x]) * 1.05
  stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
  pdf(paste0("outputs/plots_tables_objects/",patterns[x],"_cogaps.pdf"), width = 4, height = 6)
  print(ggplot(colDataDf, aes_string(color = "tissue_type", x = "tissue_type", y = patterns[x])) +
          geom_boxplot() + 
          geom_jitter(shape=16, position=position_jitter(0.2)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))+
          scale_color_manual(name = "Tissue type", labels = c("ND", "PDAC", "PNI", "VI"), values = c("#D6A857", "#DF2E2E", "#E475FF", "#2196F3")) +
          labs(color = "Tissue Type", x = "Tissue Type", y = "Pattern weight", title = paste0("Pattern ", x)) +
          theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1),
                plot.title = element_text(hjust = 0.5, size = 18)) +
          stat.test)
  dev.off()
}
#####

# GSEA and lollipop plots
# reading in hallmarks genes and cogaps patterns
all_gene_sets = msigdbr(species = "human", category = "H")
msigdbr_list = split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name)
mydata <- as.data.frame(jp4thresh$PatternMarkerScores)

# changing the sign of the pattern marker scores for illustration
mydata <- mydata * -1

# fgsea code
# performing fgsea for each CoGAPS pattern list, rank ordered by the score
mylist <- list()
for(x in 1:length(colnames(mydata))){
  tmp <- mydata[,x]
  names(tmp) <- rownames(mydata)
  tmpfgsea <- fgsea(pathways = msigdbr_list, stats = tmp)
  tmpfgsea <- as.data.frame(tmpfgsea)
  NEStmp <- tmpfgsea[order(tmpfgsea$NES, decreasing = T),]
  mylist[[x]] <- NEStmp
}

# FIGURE 3E, 4A, 5A ---- make cogaps lollipop plots, plotting the NES of each pathway ####
mylist2 <- list()
for(x in 1:length(mylist)){
  positives <- mylist[[x]] %>% dplyr::filter(NES > 0)
  negatives <- mylist[[x]] %>% dplyr::filter(NES < 0)
  positives <- positives[order(positives$NES, decreasing = T), ]
  negatives <- negatives[order(negatives$NES, decreasing = F), ]
  positives <- positives %>% dplyr::mutate(group = "positives")
  negatives <- negatives %>% dplyr::mutate(group = "negatives")
  positives <- positives[nrow(positives):1, ] 
  sig_pos <- positives %>% filter(padj < 0.05) %>% dplyr::select("pathway")
  sig_neg <- negatives %>% filter(padj < 0.05) %>% dplyr::select("pathway")
  thedata <- rbind(negatives, positives)
  #printing leading edge genes as a separate table because it is a list, which breaks the write.csv function
  theleadingedge <- t(plyr::ldply(thedata$leadingEdge, rbind))
  colnames(theleadingedge) <- thedata$pathway[!is.na(thedata$pathway)]
  thedata <- thedata %>% dplyr::select("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "group")
  thedata$pathway <- factor(thedata$pathway, levels = thedata$pathway)
  thedatatrimmed <- thedata %>% filter(NES >= 2 | NES <= -2)
  write.csv(thedata, paste0("outputs/plots_tables_objects/fgsea_results_Pattern_",x,".csv"))
  write.csv(theleadingedge, paste0("outputs/plots_tables_objects/fgsea_results_Pattern_",x,"_leading_edge.csv"))
  pdf(paste0("outputs/plots_tables_objects/lollipop_Pattern_",x,".pdf"))
  print(ggplot(thedata, aes(x=pathway, y=NES)) +
          geom_segment(aes(x=pathway, xend=pathway, y=0, yend=NES), size=0.5, alpha=0.5) +
          geom_point(size = 3, aes(fill= group), alpha=0.7, shape=21, stroke=1.4) +
          scale_fill_manual(values = c("blue", "red")) +
          theme_light() +
          theme(
            legend.position = "none",
            panel.border = element_blank(),
          ) +
          xlab("") +
          ylab("Normalized Enrichment Score") +
          coord_flip() +
          theme(axis.text.y = element_text(size = 6.8)) +
          scale_colour_manual("red", "blue"))
  dev.off()
  pdf(paste0("outputs/plots_tables_objects/trimmed_lollipop_Pattern_",x,".pdf"), height = nrow(thedatatrimmed)/4)
  print(ggplot(thedatatrimmed, aes(x=pathway, y=NES)) +
          geom_segment(aes(x=pathway, xend=pathway, y=0, yend=NES), size=0.1, alpha=0.5) +
          geom_point(size = 3, aes(fill= group), alpha=0.7, shape=21, stroke=1.4) +
          scale_fill_manual(values = c("blue", "red")) +
          theme_light() +
          theme(
            legend.position = "none",
            panel.border = element_blank(),
          ) +
          xlab("") +
          ylab("Normalized Enrichment Score") +
          coord_flip() +
          theme(axis.text.y = element_text(size = 6.8)) +
          scale_colour_manual("red", "blue"))
  dev.off()
}
#####

# VENN DIAGRAMS
# reading in the 4 CoGAPS patterns
pdac <- jp4thresh$PatternMarkers[[1]]
vi <- jp4thresh$PatternMarkers[[3]]
normal <- jp4thresh$PatternMarkers[[4]]

# reading in the VI vs PDAC DEGs
zedegs <- read.csv('outputs/plots_tables_objects/Differential_Expression_VIvsPDAC_DEGS.csv')
ups <- vivspdacdegs %>% filter(Group == "Uppers")
downs <- vivspdacdegs %>% filter(Group == "Lowers")

# getting the genes only
ups <- rownames(ups)
downs <- rownames(downs)

# taking intersection between cogaps and degs
viups <- intersect(ups, vi)
pdacdowns <- intersect(downs, pdac)
normalups <- intersect(ups, normal)

# FIGURE 4B ---- venn diagram intersection between all 138 Pattern 3 markers and all 58 genes significantly upregulated in VI relative to PDAC ####
venn.diagram(x = list(vi, ups), filename = "outputs/plots_tables_objects/vivenn.png", category.names = c("", ""),
             fill = c("#2196F3", "#2196F3"), main = "VI Enriched")
#####
# FIGURE 5C ---- venn diagram intersection between all 121 pattern 1 marker genes and all 80 genes which are significantly elevated in PDAC relative to VI ####
venn.diagram(x = list(pdac, downs), filename = "outputs/plots_tables_objects/PDACvenn.png", category.names = c("", ""),
             fill = c("#DF2E2E", "#DF2E2E"), main = "PDAC Enriched")
#####
# FIGURE 3D ---- venn diagram intersection between all 167 pattern 4 marker genes and all 58 genes which are significantly elevated in VI relative to PDAC ####
venn.diagram(x = list(normal, ups), filename = "outputs/plots_tables_objects/normalvenn.png", category.names = c("", ""),
             fill = c("#D6A857", "#2196F3"), main = "ND Enriched")
#####

# heatmap of the PDAC DEGs and VI DEGs intersected with the CoGAPS patterns
# prepping data for heatmap
significant_degs_table <- batch_adjusted[rownames(batch_adjusted) %in% c(viups, pdacdowns),]
datScale <- t(apply(significant_degs_table,1,scale))
colnames(datScale) <- colnames(batch_adjusted)
column_order = order(metadata$Type)

vivspdacdegs_1 <- vivspdacdegs[c(viups,pdacdowns),]
# making row annotation for the heatmap
groupindex <- match(rownames(datScale), c(viups, pdacdowns))
groupannonames <- vivspdacdegs_1$Group[groupindex]
row_ha = rowAnnotation(DEGs = groupannonames, col = list(DEGs = c("Uppers" = "#619CFF", "Lowers" = "#F8766D")),
                       annotation_legend_param = list(DEGs = 
                                                        list(title = "",
                                                             at = c("Lowers", "Uppers"),
                                                             labels = c("Pattern 1 markers elevated in PDAC", "Pattern 3 markers elevated in VI"))))

# FIGURE S3A ---- Heatmap of genes which are both Pattern 3 markers and significantly elevated in VI relative to PDAC ####
pdf(file = "outputs/plots_tables_objects/VIPDAC_DEGS_Patternmarkers_heatmap.pdf")
Heatmap(datScale, show_row_names = T, show_column_names = F,
        name = "Matrix",
        clustering_distance_rows = 'pearson',
        column_order = column_order,
        row_names_gp = grid::gpar(fontsize = 8, fontface = "italic"),
        top_annotation = HeatmapAnnotation(Type =
                                             metadata[,c('Type')],
                                           Patient.Alias = metadata[,c('Patient.Alias')], col = list(Type= c("ND" = "#D6A857", "PDAC" = "#DF2E2E", "PNI" = "#E475FF", "VI" = "#2196F3"),
                                                                                               Patient.Alias = c("1" = "#00BCD4", "2" = "#000000", "3" = "#CBCBCB", "4" = "#9F5500", "5" = "#7AFF88", "6" = "#3F51B5", "7" = "#FF6600", "8" = "#FFBAE2")),
                                           annotation_legend_param = list(Patient.Alias = list(title = "Patient alias"), Type = list(title = "Tissue type"))), 
        right_annotation = row_ha)
dev.off()
#####


# boxplots of individual cogaps pattern genes
# thresholded pattern markers
jp4thresh <- readRDS("outputs/thresholded_patterns.rds")
thedata <- read.csv("outputs/preprocessing/fully_batch_corrected_vsd.csv", row.names = 1)
metadata <- read.csv("outputs/preprocessing/metadata.csv", row.names = 1)
p_table <- readRDS("outputs/plots_tables_objects/p_table.rds")

# FIGURE 5D ---- Boxplots of normalized expression of genes which are both pattern 1 markers and significantly upregulated in PDAC compared to VI ####
genesSig <- jp4thresh$PatternMarkers[[1]]
metadataDat <- data.frame(metadata,t(thedata[genesSig,]))
genesSig <- genesSig[which(genesSig %in% colnames(metadataDat))]

for (g in genesSig) {
  metadataDatgene <- metadataDat %>% dplyr::select("Type", g)
  data_for_p_val_manual <- p_table %>% filter(Genes == make.names(g)) %>% dplyr::select("group_1", "group_2", "p")
  colnames(data_for_p_val_manual) <- c("group1", "group2", "p.adj")
  data_for_p_val_manual <- data_for_p_val_manual[c(1,3),]
  datman <- tibble(data_for_p_val_manual)
  
  # Convert p-values to asterisks
  datman$p.adj.signif <- sapply(datman$p.adj, p_to_asterisk)
  
  yposition <- max(metadataDatgene[,make.names(g)]) * 1.05
  stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.adj.signif", hide.ns = F, size = 6)
  pdf(paste0("outputs/geneplots/type/Pattern_1_stats/", g,"_boxplot.pdf"), width = 4, height = 6)
  print(ggplot(metadataDatgene, aes_string(x='Type', y=make.names(g), color="Type")) + 
          geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
          labs(color = "Tissue type", x = "", y = "Normalized expression", title = g) +
          theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1),
                plot.title = element_text(hjust = 0.5, size = 18, face = "italic")) +
          scale_color_manual(name = "Tissue type", labels = c("ND", "PDAC", "PNI", "VI"), values = c("#D6A857", "#DF2E2E", "#E475FF", "#2196F3")) +
          stat.test)
  dev.off()
}
#####

# FIGURE 4C ---- Boxplots of normalized expression of genes which are both pattern 3 markers and significantly upregulated in PDAC compared to VI ####
genesSig <- jp4thresh$PatternMarkers[[3]]
metadataDat <- data.frame(metadata,t(thedata[genesSig,]))
genesSig <- genesSig[which(genesSig %in% colnames(metadataDat))]
for (g in genesSig) {
  metadataDatgene <- metadataDat %>% dplyr::select("Type", g)
  data_for_p_val_manual <- p_table %>% filter(Genes == make.names(g)) %>% dplyr::select("group_1", "group_2", "p")
  colnames(data_for_p_val_manual) <- c("group1", "group2", "p.adj")
  data_for_p_val_manual <- data_for_p_val_manual[c(1,3),]
  datman <- tibble(data_for_p_val_manual)
  
  # Convert p-values to asterisks
  datman$p.adj.signif <- sapply(datman$p.adj, p_to_asterisk)
  
  datman$p.adj <- signif(datman$p.adj,digits=3)
  yposition <- max(metadataDatgene[,make.names(g)])*1.05
  stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.adj.signif", hide.ns = F, size = 6)
  pdf(paste0("outputs/geneplots/type/Pattern_3_stats/", g,"_boxplot.pdf"), width = 4, height = 6)
  print(ggplot(metadataDatgene, aes_string(x='Type', y=make.names(g), color="Type")) + 
          geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
          labs(color = "Tissue type", x = "", y = "Normalized expression", title = g) +
          theme(axis.text.x = element_text(color = "black",
                                           size = 12, angle = 45, hjust = 1),
                plot.title = element_text(hjust = 0.5, size = 18, face = "italic")) + 
          scale_color_manual(name = "Tissue type", labels = c("ND", "PDAC", "PNI", "VI"), values = c("#D6A857", "#DF2E2E", "#E475FF", "#2196F3")) +
          stat.test)
  dev.off()
}
#####
# FIGURE 3F, S2A ---- Boxplots of normalized expression of genes which are both pattern 4 markers and significantly upregulated in PDAC compared to VI ####
genesSig <- jp4thresh$PatternMarkers[[4]]
metadataDat <- data.frame(metadata,t(thedata[genesSig,]))
genesSig <- genesSig[which(genesSig %in% colnames(metadataDat))]
for (g in genesSig) {
  metadataDatgene <- metadataDat %>% dplyr::select("Type", g)
  data_for_p_val_manual <- p_table %>% filter(Genes == make.names(g)) %>% dplyr::select("group_1", "group_2", "p")
  colnames(data_for_p_val_manual) <- c("group1", "group2", "p.adj")
  data_for_p_val_manual <- data_for_p_val_manual[c(1,3),]
  datman <- tibble(data_for_p_val_manual)
  
  # Convert p-values to asterisks
  datman$p.adj.signif <- sapply(datman$p.adj, p_to_asterisk)
  
  datman$p.adj <- signif(datman$p.adj,digits=3)
  yposition <- max(metadataDatgene[,make.names(g)]) * 1.05
  stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.adj.signif", hide.ns = F, size = 6)
  pdf(paste0("outputs/geneplots/type/Pattern_4_stats/", g,"_boxplot.pdf"), width = 4, height = 6)
  print(ggplot(metadataDatgene, aes_string(x='Type', y=make.names(g), color="Type")) + 
          geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
          labs(color = "Tissue type", x = "", y = "Normalized expression", title = g) +
          theme(axis.text.x = element_text(color = "black",
                                           size = 12, angle = 45, hjust = 1),
                plot.title = element_text(hjust = 0.5, size = 18, face = "italic")) + 
          scale_color_manual(name = "Tissue type", labels = c("ND", "PDAC", "PNI", "VI"), values = c("#D6A857", "#DF2E2E", "#E475FF", "#2196F3")) +
          stat.test)
  dev.off()
}
#####

# COGAPS BARCODE PLOTS
# HALLMARK EMT COGAPS 1 BARCODE PLOT
# reading in the HALLMARK gene sets
all_gene_sets = msigdbr(species = "human", category = "H")
msigdbr_list = split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name)
mydata <- as.data.frame(jp4thresh$PatternMarkerScores)
mydata <- mydata * -1
emt <- msigdbr_list$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
index_emt1 <- as.vector(na.omit(match(emt, rownames(mydata))))
# FIGURE 5B ---- Barcode plot of the genes within the HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION ####
pdf("outputs/plots_tables_objects/emt_barcode_pattern_1.pdf")
# setting quantiles to be extremely large so that the barcode plot isn't arbitrarily colored blue/grey/red. 
barcodeplot(mydata$Pattern_1, index_emt1, labels = c("Low", "High"), xlab = "Pattern 1 Marker Score", quantiles = c(-9999,9999)*sqrt(9999))
dev.off()
#####

# Classical and basal barcode plots
# reading in all the genes in VIxPDAC ranked by log fold change
vipdac_all <- read.csv("outputs/plots_tables_objects/Differential_Expression_VIvsPDAC_gsea_allgenes_rankedLFC.csv")

# STORING THE GENE SETS FOR CLASSICAL, BASAL, MESENCHYMAL, AND COLLECTIVE
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

index_collective <- as.vector(na.omit(match(collective, vipdac_all$rn)))
index_mesenchymal <- as.vector(na.omit(match(mesenchymal, vipdac_all$rn)))
newstat <- vipdac_all$logFC

# statistics for each of the barcode plots
geneSetTest(index_collective, newstat, alternative = "down", ranks.only = F)
p_mesenchymal <-  geneSetTest(index_mesenchymal, newstat, alternative = "down", ranks.only = F)
p_collective <- geneSetTest(index_collective, newstat, alternative = "up", ranks.only = F)
geneSetTest(index_mesenchymal, newstat, alternative = "up", ranks.only = F)

index_classical <- as.vector(na.omit(match(classical, vipdac_all$rn)))
index_basal <- as.vector(na.omit(match(basal, vipdac_all$rn)))
newstat <- vipdac_all$logFC

# statistics for each of the barcode plots
geneSetTest(index_classical, newstat, alternative = "down", ranks.only = F)
p_basal <- geneSetTest(index_basal, newstat, alternative = "down", ranks.only = F)
p_classical <- geneSetTest(index_classical, newstat, alternative = "up", ranks.only = F)
geneSetTest(index_basal, newstat, alternative = "up", ranks.only = F)

# Making the barcode plots
# FIGURE 5E ---- Barcode plots of marker genes from the Moffitt classification phenotypes. ####
pdf("outputs/plots_tables_objects/barcodeplot_classical.pdf")
# setting quantiles to be extremely large so that the barcode plot isn't arbitrarily colored blue/grey/red. 
barcodeplot(newstat, index_classical, labels = c("PDAC", "VI"), xlab = "log2 fold change") + title(paste0("p = ", p_classical), quantiles = c(-9999,9999)*sqrt(9999))
dev.off()

pdf("outputs/plots_tables_objects/barcodeplot_basal.pdf")
# setting quantiles to be extremely large so that the barcode plot isn't arbitrarily colored blue/grey/red. 
barcodeplot(newstat, index_basal, labels = c("PDAC", "VI"), xlab = "log2 fold change") + title(paste0("p = ", p_basal), quantiles = c(-9999,9999)*sqrt(9999))
dev.off()
#####

# FIGURE S5C ---- Barcode plots of marker genes from the collective and mesenchymal invasion phenotypes. ####
pdf("outputs/plots_tables_objects/barcodeplot_collective.pdf")
# setting quantiles to be extremely large so that the barcode plot isn't arbitrarily colored blue/grey/red. 
barcodeplot(newstat, index_collective, labels = c("PDAC", "VI"), xlab = "log2 fold change") + title(paste0("p = ", p_collective), quantiles = c(-9999,9999)*sqrt(9999))
dev.off()

pdf("outputs/plots_tables_objects/barcodeplot_mesenchymal.pdf")
# setting quantiles to be extremely large so that the barcode plot isn't arbitrarily colored blue/grey/red. 
barcodeplot(newstat, index_mesenchymal, labels = c("PDAC", "VI"), xlab = "log2 fold change") + title(paste0("p = ", p_mesenchymal), quantiles = c(-9999,9999)*sqrt(9999))
dev.off()
#####

# doing GSVA for each of the 4 gene sets within each individual patient
genesets <- list(basal, classical)
names(genesets) <- c("basal", "classical")

# reading in pre batch corrected data
metadata <- read.csv("outputs/preprocessing/metadata.csv", row.names = 1)
thedata <- read.csv("outputs/preprocessing/vsd_not_batch_corrected_vsd.csv", row.names = 1)

# preparing the data for individual patient analysis of classical and basal-like
theslides <- unique(metadata$Slide)
metadataslidelist <- list()
for(x in 1:length(theslides)){
  thenumbers <- which(metadata$Slide == theslides[x] & (metadata$Type == "VI" | metadata$Type == "PDAC"))
  thedata <- as.matrix(thedata)
  thedataslide <- thedata[,thenumbers]
  metadataslide <- metadata[thenumbers,]
  theoutput <- gsva(thedataslide, genesets, verbose=FALSE)
  yeet <- theoutput[1,] > theoutput[2,]
  yeet[which(yeet == TRUE)] <- "basal"
  yeet[which(yeet == "FALSE")] <- "classical"
  metadataslide$basal_classical <- yeet
  yeetmetadataslide <- metadataslide[,c("basal_classical", "Type", "Patient.Alias")]
  metadataslide$basal_score <- theoutput[1,]
  metadataslide$classical_score <- theoutput[2,]
  metadataslidelist[[x]] <- metadataslide[,c('Name', 'basal_classical', 'basal_score', 'classical_score', 'Type', 'Slide', 'Patient.Alias')]
  
}

# importing the required library
thedata <- do.call(rbind, metadataslidelist)

# melting the data 
thedata <- reshape2::melt(thedata, id.var = c('Name', 'Type', "Slide", "Patient.Alias", 'basal_classical'),
                          variable.name = 'score')

# factoring the scores and tissue types into levels
thedata$score <- factor(thedata$score, levels = c("classical_score", "basal_score"))
thedata$Type <- factor(thedata$Type, levels = c("PDAC", "VI"))

# making p value lists for each patient comparison
classicallist <- list()
classicalnames <- list()
basallist <- list()
basalnames <- list()
for(x in 1:8){
  viclassical <- thedata %>% filter(Patient.Alias == x & Type == "VI" & score == "classical_score") %>% dplyr::select("value") %>% unlist()
  pdacclassical <- thedata %>% filter(Patient.Alias == x & Type == "PDAC" & score == "classical_score") %>% dplyr::select("value") %>% unlist()
  vibasal <- thedata %>% filter(Patient.Alias == x & Type == "VI" & score == "basal_score") %>% dplyr::select("value") %>% unlist()
  pdacbasal <- thedata %>% filter(Patient.Alias == x & Type == "PDAC" & score == "basal_score") %>% dplyr::select("value") %>% unlist()
  classical_p <- wilcox.test(viclassical, pdacclassical)
  basal_p <- wilcox.test(vibasal, pdacbasal)
  classicallist[[x]] <- classical_p$p.value
  classicalnames[[x]] <- x
  basallist[[x]] <- basal_p$p.value
  basalnames[[x]] <- x
}

theclassical <- data.frame(score = unlist(classicallist), group = unlist(classicalnames))
thebasal <- data.frame(score = unlist(basallist), group = unlist(basalnames))
thetable <- rbind(theclassical, thebasal)
thetable$scoregroup <- c(rep("classical", 8), rep('basal', 8))

# FIGURE S5A ---- Comparison of VI (blue) and PDAC (red) by gene set variation analysis (GSVA) of classical and basal-like phenotypes ####
# making figure for intra-patient comparisons of classical and basal-like scores
pdf("outputs/plots_tables_objects/classical_basal_facet.pdf", width = 6, height = 12)
ggplot(thedata, aes_string(x='score', y="value", fill='Type')) + 
  geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=16, size=2, color = 'gold',
                                position = position_dodge2(width = 0.75,   
                                                           preserve = "single")) +
  geom_jitter(shape=16, position=position_dodge(0.75)) +
  facet_wrap(~Patient.Alias, ncol =1) +
  scale_x_discrete(name = "Moffitt Classification", labels = c("Classical", "Basal-like")) +
  scale_fill_manual(name = "Tissue Type", labels = c("PDAC", "VI"), values = c("#DF2E2E", "#2196F3")) + 
  scale_y_continuous(name = "GSVA Enrichment Score") + 
  theme_bw() +
  theme(panel.background=element_rect(colour="black")) +
  geom_vline(xintercept = 1.5, color = "black") +
  #ylim(-0.5, 0.6) +
  guides(fill = guide_legend(override.aes = list(shape = c(NA, NA))))
dev.off()
#####

# making the gene set lists for mesenchymal and collective
genesets <- list(mesenchymal, collective)
names(genesets) <- c("mesenchymal", "collective")
theslides <- unique(metadata$Slide)
metadataslidelist <- list()
metadata <- read.csv("outputs/preprocessing/metadata.csv", row.names = 1)
thedata <- read.csv("outputs/preprocessing/vsd_not_batch_corrected_vsd.csv", row.names = 1)

theslides <- unique(metadata$Slide)
metadataslidelist <- list()

# splitting the data into individual patients
for(x in 1:length(theslides)){
  thenumbers <- which(metadata$Slide == theslides[x] & (metadata$Type == "VI" | metadata$Type == "PDAC"))
  thedatamat <- as.matrix(thedata)
  thedataslide <- thedatamat[,thenumbers]
  metadataslide <- metadata[thenumbers,]
  theoutput <- gsva(thedataslide, genesets, verbose=FALSE)
  yeet <- theoutput[1,] > theoutput[2,]
  yeet[which(yeet == TRUE)] <- "collective"
  yeet[which(yeet == "FALSE")] <- "mesenchymal"
  metadataslide$mesenchymal_collective <- yeet
  yeetmetadataslide <- metadataslide[,c("mesenchymal_collective", "Type", "Patient.Alias")]
  metadataslide$mesenchymal_score <- theoutput[1,]
  metadataslide$collective_score <- theoutput[2,]
  metadataslidelist[[x]] <- metadataslide[,c('Name', 'mesenchymal_collective', 'mesenchymal_score', 'collective_score', 'Type', 'Slide', 'Patient.Alias')]
  
}

thedata <- do.call(rbind, metadataslidelist)

thedata <- reshape2::melt(thedata, id.var = c('Name', 'Type', "Slide", "Patient.Alias", 'mesenchymal_collective'),
                          variable.name = 'score')


thedata$score <- factor(thedata$score, levels = c("collective_score", "mesenchymal_score"))
thedata$Type <- factor(thedata$Type, levels = c("PDAC", "VI"))

# making p values for each comparison
collectivelist <- list()
collectivenames <- list()
mesenchymallist <- list()
mesenchymalnames <- list()
for(x in 1:8){
  vicollective <- thedata %>% filter(Patient.Alias == x & Type == "VI" & score == "collective_score") %>% dplyr::select("value") %>% unlist()
  pdaccollective <- thedata %>% filter(Patient.Alias == x & Type == "PDAC" & score == "collective_score") %>% dplyr::select("value") %>% unlist()
  vimesenchymal <- thedata %>% filter(Patient.Alias == x & Type == "VI" & score == "mesenchymal_score") %>% dplyr::select("value") %>% unlist()
  pdacmesenchymal <- thedata %>% filter(Patient.Alias == x & Type == "PDAC" & score == "mesenchymal_score") %>% dplyr::select("value") %>% unlist()
  collective_p <- wilcox.test(vicollective, pdaccollective)
  mesenchymal_p <- wilcox.test(vimesenchymal, pdacmesenchymal)
  collectivelist[[x]] <- collective_p$p.value
  collectivenames[[x]] <- x
  mesenchymallist[[x]] <- mesenchymal_p$p.value
  mesenchymalnames[[x]] <- x
}


thecollective <- data.frame(score = unlist(collectivelist), group = unlist(collectivenames))
themesenchymal <- data.frame(score = unlist(mesenchymallist), group = unlist(mesenchymalnames))
thetable <- rbind(thecollective, themesenchymal)
thetable$scoregroup <- c(rep("collective", 8), rep('mesenchymal', 8))

# FIGURE S5B collective and mesenchymal invasion phenotypes within individual patients ####
pdf("outputs/plots_tables_objects/collective_mesenchymal_facet.pdf", width = 6, height = 12)
ggplot(thedata, aes_string(x='score', y="value", fill='Type')) + 
  geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=16, size=2, color = 'gold',
                                position = position_dodge2(width = 0.75,   
                                                           preserve = "single")) +
  geom_jitter(shape=16, position=position_dodge(0.75)) +
  facet_wrap(~Patient.Alias, ncol =1) +
  scale_x_discrete(name = "Invasion Phenotype", labels = c("Collective", "Mesenchymal")) +
  scale_fill_manual(name = "Tissue Type", labels = c("PDAC", "VI"), values = c("#DF2E2E", "#2196F3")) + 
  scale_y_continuous(name = "GSVA Enrichment Score") + 
  theme_bw() +
  theme(panel.background=element_rect(colour="black")) +
  geom_vline(xintercept = 1.5, color = "black") +
  guides(fill = guide_legend(override.aes = list(shape = c(NA, NA))))
dev.off()
#####
