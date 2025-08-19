library(ggplot2)
library(ggpubr)
library(CoGAPS)
library(plyr)
library(dplyr)
library(forcats)
library(monocle3)
library(Seurat)
library(projectR)
library(dittoSeq)
library(MAST)
library(hdf5r)

source("scripts/00_Custom_Functions_atlas.R")

# load in PDAC atlas
thedata <- readRDS("inputs/cds_combined_epithelial.rds")

# load in the GeoMx CoGAPS object
cogaps <- readRDS("inputs/batch_merged_cogapsP4.rds")

# Assign cells as pattern markers based on the patternMarkers function
## Run the pattern marker statistic on the cells 
cellPatternMarker <- patternMarkers(cogaps, axis = 2)

## Add the pattern to the cell marker result
for(i in 1:length(cellPatternMarker$PatternMarkers)){
  ## Concert to a df
  cellPatternMarker$PatternMarkers[[i]] <- as.data.frame(cellPatternMarker$PatternMarkers[[i]])
  ## create column name
  colnames(cellPatternMarker$PatternMarkers[[i]]) <- "CellID"
  cellPatternMarker$PatternMarkers[[i]]$Pattern <- names(cellPatternMarker$PatternMarkers[i])
}

## Combined the tables
cellPatternMarkerDataFrame <- do.call("rbind", cellPatternMarker$PatternMarkers)
## Update rownames
rownames(cellPatternMarkerDataFrame) <- cellPatternMarkerDataFrame$CellID

## Add the pattern assignment to the thedata
colData(thedata)$AssignedPattern <- cellPatternMarkerDataFrame[colnames(thedata), "Pattern"]

# create an additional classifier for cells that divides Epithelial_normal
# into normal cells from non-PDAC patients and normal_tumor_adjacent cells
# from patients with PDAC
thedata@colData$Epithelial_Cell_Type <- paste0(as.character(thedata@colData$TN),"_", 
                                               as.character(thedata@colData$TN_assigned_cell_type_immune_broad))
thedata@colData$Epithelial_Cell_Type <- gsub("T_Epithelial_cancer", "Cancer",
                                             gsub("T_Epithelial_unspecified", "Unspecified",
                                                  gsub("T_Epithelial_normal", "Normal_Tumor_Adjacent",
                                                       gsub("N_Epithelial_cancer", "Cancer",
                                                            gsub("N_Epithelial_unspecified", "Unspecified",
                                                                 gsub("N_Epithelial_normal", "Normal", thedata@colData$Epithelial_Cell_Type))))))
# change row names from ensembl to hugo symbol
rownames(thedata) <- fData(thedata)$gene_short_name

# S method for SingleCellExperiment
thedata <- as.Seurat(
  thedata,
  counts = "counts",
  data = NULL,
  assay = NULL,
  project = "SingleCellExperiment"
)

# normalizing the data
thedata <- NormalizeData(
  thedata,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1,
  block.size = NULL,
  verbose = TRUE,
)

# this step regresses out the manuscript source (eg peng, steele) and scales the data
# commented out because this step takes absolutely forever, and I've already done it once, so I just read it in again below
#thedata <- ScaleData(
#   thedata,
#   assay = "originalexp",
#   vars.to.regress = "manuscript"
# )

#saveRDS(thedata@assays$originalexp@scale.data, "outputs/plots_tables_objects/pdac_atlas_scaledata.rds")
thedata@assays$originalexp@scale.data <- readRDS("outputs/plots_tables_objects/pdac_atlas_scaledata.rds")

mydf <- thedata@assays$originalexp@scale.data

# project CoGAPS pattern weights from PDAC atlas epithelial cells onto
# SCT scaled transcript matrix of panin visium epithelial cells
project_embed <- projectR(data = as.matrix(mydf), loadings = cogaps)

# add projected patterns to seurat object
pattern_weights <- as.data.frame(t(project_embed))

# list of CoGAPS Patterns
patterns <- colnames(pattern_weights)

# adding the patterns to the object metadata
for(pattern in colnames(pattern_weights)) {
  thedata <- AddMetaData(object = thedata,
                         metadata = pattern_weights[[pattern]],
                         col.name = pattern)
}

# getting the patterns from the object metadata
thepatterns <- thedata@meta.data[, 27:30]

#making a loop that labels cells based on the pattern that they express most intensely
the_vector <- c()
for(x in 1:nrow(thepatterns)){
  maxval <- max(thepatterns[x,])
  max_index <- which(thepatterns[x,] == maxval)
  the_vector[x] <- colnames(thepatterns)[max_index]
}

# storing the most intensely expressed pattern
thedata$patterns <- the_vector

# getting the cancerous cells
cancers <- thedata@meta.data %>% filter(TN_assigned_cell_type_immune_broad != "Epithelial_normal")

# getting the third quartile score for each pattern
sum1 <- summary(cancers$Pattern_1)
sum3 <- summary(cancers$Pattern_3)
sum4 <- summary(cancers$Pattern_4)
q3_1 <- sum1[5]
q3_3 <- sum3[5]
q3_4 <- sum4[5]

# getting cells which score above the third quartile for each pattern and score higher than the other two patterns
thedata$geomx_labels <- "unlabeled"
PDAC_index <- which(thedata@meta.data$Pattern_1 > q3_1 & thedata@meta.data$patterns == "Pattern_1")
VI_index <- which(thedata@meta.data$Pattern_3 > q3_3 & thedata@meta.data$patterns == "Pattern_3")
Normal_index <- which(thedata@meta.data$Pattern_4 > q3_4 & thedata@meta.data$patterns == "Pattern_4")

# labeling the cells which score above the third quartile for patterns 1 and 3 as Projected PDAC and Projected VI, respectively 
thedata@meta.data$geomx_labels <- "Undefined"
thedata@meta.data$geomx_labels[PDAC_index] <- "Projected PDAC"
thedata@meta.data$geomx_labels[VI_index] <- "Projected VI"
thedata@meta.data$geomx_labels[Normal_index] <- "Projected ND"

# factoring the labels
thedata@meta.data$geomx_labels <- factor(thedata@meta.data$geomx_labels, levels=c('Projected ND', "Projected PDAC", "Projected VI", "Undefined"))

# MAKING THE VIOLIN PLOTS OF THE OVERLAPPING GENES BETWEEN GEOMX DE AND PROJECTR ATLAS DE
#markers <- FindMarkers(thedata, group.by = "geomx_labels", ident.1 = "Projected VI", ident.2 = "Projected PDAC", test.use = "MAST", slot = "scale")

#saveRDS(markers, "outputs/plots_tables_objects/PengSteeleMarkers.rds")

themarkers <- readRDS("outputs/plots_tables_objects/PengSteeleMarkers.rds")
vi.vs.pdac <- readRDS("outputs/plots_tables_objects/vipdaclimma.rds")

# #filtering out significant genes (used 0.58 as lfc cutoff for volcano plots)
vi.vs.pdac %>% filter(adj.P.Val < 0.05 & (logFC > 0.58 | logFC < -0.58)) -> significant_degs

# replacing zero with lowest number R can represent
themarkers[which(themarkers$p_val_adj == 0), 5] <- .Machine$double.xmin

# getting significant degs
significant_atlas <- themarkers %>% filter(p_val_adj < 0.05 & (avg_log2FC > 0.58 | avg_log2FC < -0.58))

# Getting positives and negatives for GeoMx DEGs and Atlas DEGs
geomx_neg <- significant_degs %>% filter(logFC < 0)
geomx_pos <- significant_degs %>% filter(logFC > 0)

atlas_neg <- significant_atlas %>% filter(avg_log2FC < 0)
atlas_pos <- significant_atlas %>% filter(avg_log2FC > 0)

# Getting the intersection between the two groups
negatives <- intersect(rownames(atlas_neg), rownames(geomx_neg))
positives <- intersect(rownames(atlas_pos), rownames(geomx_pos))

themarkers$group <- "unchanged"
themarkers$group[rownames(themarkers) %in% rownames(atlas_neg)] <- "down"
themarkers$group[rownames(themarkers) %in% rownames(atlas_pos)] <- "up"
write.csv(themarkers, "outputs/plots_tables_objects/Projected_VI_vs_Projected_PDAC_results.csv")
thedata <- SetIdent(thedata, value = "geomx_labels")

# MAKING HEATMAP OF THE POSITIVES AND NEGATIVES
thegenes <- c(positives, negatives)
metadata1 <- thedata@meta.data
thecells <- as.vector(unlist(CellsByIdentities(thedata, idents = c("Projected PDAC", "Projected VI"), cells = NULL)))

# FIGURE S4C ####
# making heatmap of the projected VI and PDAC populations and the genes which overlap between the DEGs from the projected populations and the DEGs from the GeoMx comparison
pdf("outputs/plots_tables_objects/atlas_heatmap_overlap_of_geomx_and_seurat.pdf")
DoHeatmap(thedata, features = thegenes, group.by = "geomx_labels", cells = thecells, label = F, group.colors = c("#D6A857", "#DF2E2E", "#2196F3", "#CBCBCB")) + 
  theme(text = element_text(size = 8),
        axis.text.y = element_text(face = "italic"))
dev.off()
#####

thedata@meta.data$TN_assigned_cell_type_immune_specific <- factor(thedata@meta.data$TN_assigned_cell_type_immune_specific, levels = c("Epithelial_normal", "Epithelial_cancer", "Epithelial_unspecified"))

# FIGURE S4B ####
pdf("outputs/plots_tables_objects/umap_projected_cell_types.pdf")
DimPlot(thedata, group.by = "geomx_labels", cols = c("#D6A857", "#DF2E2E", "#2196F3", "#CBCBCB"))
dev.off()
#####

# FIGURE S4A ####
# take pictures of pattern umaps
pdf("outputs/plots_tables_objects/atlas_projection_pattern_1.pdf")
FeaturePlot(thedata, features = "Pattern_1", slot = "scale")
dev.off()
pdf("outputs/plots_tables_objects/atlas_projection_pattern_2.pdf")
FeaturePlot(thedata, features = "Pattern_2", slot = "scale")
dev.off()
pdf("outputs/plots_tables_objects/atlas_projection_pattern_3.pdf")
FeaturePlot(thedata, features = "Pattern_3", slot = "scale")
dev.off()
pdf("outputs/plots_tables_objects/atlas_projection_pattern_4.pdf")
FeaturePlot(thedata, features = "Pattern_4", slot = "scale")
dev.off()
#####
