# Loading packages
library(tidyverse)
library(dplyr)
library(epitools)
library(ggpubr)
library(ggprism)
library(shiny)
library(ggplot2)
sessionInfo()

source("scripts/00_Custom_Functions.R")

# Importing HALO outputs
combined_table <- read.csv("inputs/raw_HALO_outputs.csv", row.names = 1)
unique(combined_table$Image.Location)
# Filtering for epithelial cells
combined_table_epi <- combined_table %>% filter(PanCK.Positive.Classification == 1)

# Selecting epithelial cells for being within or outside of the PDAC_enriched annotations
the_pdac <- combined_table_epi %>% filter(Analysis.Region == "PDAC_enriched")
the_other <- combined_table_epi %>% filter(Analysis.Region != "PDAC_enriched")

# shortening the table for computational efficiency of next step
the_other_short <- the_other %>% dplyr::select("Image.Location", "XMin", "XMax", "YMin", "YMax")
the_pdac_short <- the_pdac %>% dplyr::select("Image.Location", "XMin", "XMax", "YMin", "YMax")

# Because some of the "other" (ie PNI, ND, duodenum) regions are contained within PDAC enriched regions, they are double counted. So I am anti-joining the PDAC enriched with the non-PDAC enriched area (return all PDAC enriched rows without matches in the non-PDAC enriched)
the_pdac_short_distinct <- anti_join(the_pdac_short, the_other_short)

# By joining the unique PDAC cells with all of the other cells, I get a complete table of all epithelial cells in the data without any duplicates
deduplicated_cells <- rbind(the_pdac_short_distinct, the_other_short)

# getting the table with all the data subset on the non-duplicated cells
deduplicated_table <- combined_table_epi[rownames(deduplicated_cells),]

# Filtering out normal duct and VI conventional for downstream comparisons
deduplicated_table <- deduplicated_table %>% filter(Analysis.Region %in% c("PDAC_enriched", "VI_IN_LIKE", "VI_DESTRUCTIVE"))
test_for_duplicated <- deduplicated_table %>% dplyr::select("Image.Location", "XMin", "XMax", "YMin", "YMax")

# Reading in the genes that will be used in downstream analysis
the_markers <- c("DAPI", "TFF1", "CLDN18.2", "GATA6", "KRT5", "KRT17", "S100A2", "CEACAM5", "MUC13", "LAMC2")

# testing for duplicates
which(duplicated(test_for_duplicated))
colnames(deduplicated_table)
# Cutting the individual slides into equal sized tiles (500x500 pixel), getting cells only in the PDAC_enriched areas
output_by_slide <- list()
final_table_list <- list()
for(z in 1:length(unique(deduplicated_table$Image.Location))){
  theimage <- unique(deduplicated_table$Image.Location)[z]
  theindex <- which(deduplicated_table$Image.Location == theimage)
  tmp_table <- deduplicated_table[theindex,]
  tmp_table <- tmp_table %>% filter(PanCK.Positive.Classification == 1 & Analysis.Region == "PDAC_enriched")
  
  xmin <- tmp_table %>% dplyr::select("XMin") %>% min() 
  xmax <- tmp_table %>% dplyr::select("XMax") %>% max()
  ymin <- tmp_table %>% dplyr::select("YMin") %>% min()
  ymax <- tmp_table %>% dplyr::select("YMax") %>% max()
  
  mincoord <- c(xmin, ymin)
  maxcoord <- c(xmax, ymax)
  
  xlen <- maxcoord[1] - mincoord[1]
  ylen <- maxcoord[2] - mincoord[2]
  tile_xdim <- 500
  tile_ydim <- 500
  
  the_tile_lists_by_x <- list()
  the_tile_tables_by_y <- list()
  thecount = 0
  xloop <- ceiling(xlen/tile_xdim)
  yloop <- ceiling(ylen/tile_ydim)
  # Looping through all tile coordinates. The four elements of the vector represent x1,y1,x2,y2
  for(x in 1:xloop){
    the_tile_tables_by_y <- list()
    for(y in 1:yloop){
      thecount <- thecount+1
      # This makes the tile coordinates. The four elements of the vector represent x1,y1,x2,y2
      tile_coords <- c(mincoord[1] + tile_xdim*(x-1), mincoord[2] + tile_ydim*(y-1), mincoord[1] + tile_xdim*x, mincoord[2] + tile_ydim*y)
      
      # getting the centroid coordinates for each cell
      tmp_table$XMean <- tmp_table %>% dplyr::select("XMin", "XMax") %>% rowMeans()
      tmp_table$YMean <- tmp_table %>% dplyr::select("YMin", "YMax") %>% rowMeans()
      
      # this gets all of the cells within the tile
      tmp <- tmp_table %>% filter((XMean > tile_coords[1] & XMean < tile_coords[3]) & (YMean > tile_coords[2] & YMean < tile_coords[4]))
      
      #only using tiles if there are at least 10 cells in that tile
      if(nrow(tmp)  > 9) {
        tmp$tile_number <- thecount
        the_tile_tables_by_y[[y]] <- tmp
        names(the_tile_tables_by_y)[y] <- thecount
      }
      else{
        next
      }
    }
    the_tile_lists_by_x[[x]] <- do.call(rbind, the_tile_tables_by_y)
  }
  print("finished")
  final_table_list[[z]] <- do.call(rbind, the_tile_lists_by_x)
}
colnames(deduplicated_table)
# combining the processed tables with the tiles
# There are 11596 cells that are lost because they are part of a tile that has fewer than 10 cells
theresult <- do.call(rbind, final_table_list)

# Creating a column for tile-patient
theresult$tile_number_id <- paste0(theresult$tile_number,"_",theresult$Image.Location)

# Shortening the sample names
new_rownames <- sub(".*VI", "VI", rownames(theresult))

# Creating a common column between the old table (deduplicated_table, no tiles annotated) and new table (theresult, with annotated tiles)
theresult$thecellnames <- new_rownames
deduplicated_table$thecellnames <- rownames(deduplicated_table)

# Creaing notin function
`%notin%` <- negate(`%in%`)

# getting tables for each VI subtype. These tables were generated through the rshiny function where I manually named the individual VI foci
codex_01 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_01_groups.rds")
codex_05 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_05_groups.rds")
codex_09 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_09_groups.rds")
codex_10 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_10_groups.rds")
codex_12 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_12_groups.rds")
codex_14 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_14_groups.rds")
codex_16 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_16_groups.rds")
codex_17 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_17_groups.rds")
codex_19 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_19_groups.rds")
codex_21 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_21_groups.rds")
codex_23 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_23_groups.rds")
codex_25 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_25_groups.rds")
codex_26 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_26_groups.rds")
codex_28 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_28_groups.rds")
codex_29 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_29_groups.rds")
codex_32 <- readRDS("outputs/plots_tables_objects/phenocycler/codex_32_groups.rds")
KH_002 <- readRDS("outputs/plots_tables_objects/phenocycler/KH_002_groups.rds")
KH_003 <- readRDS("outputs/plots_tables_objects/phenocycler/KH_003_groups.rds")
KH_005 <- readRDS("outputs/plots_tables_objects/phenocycler/KH_005_groups.rds")

# Making a list of the named individual VI foci
groupslist <- list(codex_01, codex_05, codex_09, codex_10, codex_12,
                   codex_14, codex_16, codex_17, codex_19, codex_21, 
                   codex_23, codex_25, codex_26, codex_28, codex_29, 
                   codex_32, KH_002, KH_003, KH_005)

# Making a list that contains one patient per list element
rbind_list <- list()
for(x in 1:length(groupslist)){
  tmp <- do.call(rbind, groupslist[[x]])
  rbind_list[[x]] <- tmp
}

# Making table that contains all of the patients
full_vi_table <- do.call(rbind, rbind_list)

# Selecting only the VI, because I also annotated ND and other things. 
full_vi_table <- full_vi_table %>% filter(Analysis.Region %in% c("VI_IN_LIKE", "VI_DESTRUCTIVE"))

# cleaning row names
cleaned_text <- sub(".*VI_", "VI_", rownames(full_vi_table))
rownames(full_vi_table) <- cleaned_text
result <- sub("(\\.qptiff).*", "\\1", cleaned_text)
full_vi_table$image <- result

# Merging the deduplicated table (no PDAC foci) with the table containing the individually labeled VI foci
deduplicated_table_with_vi_foci <- merge(full_vi_table, deduplicated_table, by = "row.names")

# Removing unnecessary columns
deduplicated_table_with_vi_foci <- subset(deduplicated_table_with_vi_foci, select=-c(xmean,ymean, Analysis.Region.x, Row.names))

# Renaming the rows with the shorter row names
rownames(deduplicated_table_with_vi_foci) <- deduplicated_table_with_vi_foci$thecellnames

# Renaming columns for clarity
colnames(deduplicated_table_with_vi_foci)[which(colnames(deduplicated_table_with_vi_foci) == "name")] <- "ROI"
colnames(deduplicated_table_with_vi_foci)[which(colnames(deduplicated_table_with_vi_foci) == "Analysis.Region.y")] <- "Analysis.Region"

# Confirming that we only have VI in this table
unique(deduplicated_table_with_vi_foci$Analysis.Region)

# Making a roi-slide column (so that we can distinguish individual VI foci from different patients)
deduplicated_table_with_vi_foci$ROI_by_slide <- paste0(deduplicated_table_with_vi_foci$ROI, "_", deduplicated_table_with_vi_foci$Image.Location)

# Adding the PDAC_enriched slide ROIs/Foci to the original data, to eventually have one data frame that contains the labeled individual VI foci and PDAC foci
deduplicated_table_with_pdac_foci <- merge(theresult, deduplicated_table, by = "thecellnames")
colnames(deduplicated_table_with_pdac_foci)

# Removing duplicated/unnecessary columns created from the merge
deduplicated_table_with_pdac_foci <- deduplicated_table_with_pdac_foci[,1:66] 
deduplicated_table_with_pdac_foci <- subset(deduplicated_table_with_pdac_foci, select=-c(XMean,YMean))

# Creating matching column names
colnames(deduplicated_table_with_pdac_foci)[which(colnames(deduplicated_table_with_pdac_foci) == "tile_number")] <- "ROI"
colnames(deduplicated_table_with_pdac_foci)[which(colnames(deduplicated_table_with_pdac_foci) == "tile_number_id")] <- "ROI_by_slide"

# Matching columns between the two tables
pdac_values <- deduplicated_table_with_pdac_foci[,6:62]
vi_values <- deduplicated_table_with_vi_foci[,7:63]

# Matching metadata between the two tables
pdac_metadata <- deduplicated_table_with_pdac_foci[,c(1:5,63,64)]
vi_metadata <- deduplicated_table_with_vi_foci[,c(1:6,64,65)]
colnames(vi_metadata)[1] <- "ROI"
vi_metadata <- subset(vi_metadata, select = -image)
vi_metadata_2 <- vi_metadata
vi_metadata[,1] <- vi_metadata_2$thecellnames
colnames(vi_metadata)[1] <- "thecellnames"
vi_metadata[,6] <- vi_metadata_2$ROI
colnames(vi_metadata)[6] <- "ROI"
colnames(pdac_metadata) <- colnames(vi_metadata)
colnames(pdac_values) <- colnames(vi_values)

# Making table for annotated VI and annotated PDAC, with columns in the same order
annotated_VI <- cbind(vi_metadata, vi_values)
annotated_PDAC <- cbind(pdac_metadata, pdac_values)

# Making table that has individually named VI and individually named PDAC foci
annotated_total <- rbind(annotated_VI, annotated_PDAC)
colnames(annotated_total)[3] <- "Analysis.Region"

# All cells % positive and negative for each marker
inlike <- c()
dest <- c()
pdac <- c()

# getting the percent % for each marker by tissue type (analysis region)
for(x in 1:length(the_markers)){
  if(the_markers[x] == "CLDN18.2")
  {
    inlike[x] <- deduplicated_table %>% filter(Analysis.Region == "VI_IN_LIKE" & PanCK.Positive.Classification == 1) %>% dplyr::select(paste0(the_markers[x],".Positive.Cytoplasm.Classification")) %>% as.vector() %>% unlist() %>% as.vector() %>% mean() * 100
    dest[x] <- deduplicated_table %>% filter(Analysis.Region == "VI_DESTRUCTIVE" & PanCK.Positive.Classification == 1) %>% dplyr::select(paste0(the_markers[x],".Positive.Cytoplasm.Classification")) %>% as.vector() %>% unlist() %>% as.vector() %>% mean() * 100
    pdac[x] <- deduplicated_table %>% filter(Analysis.Region == "PDAC_enriched" & PanCK.Positive.Classification == 1) %>% dplyr::select(paste0(the_markers[x],".Positive.Cytoplasm.Classification")) %>% as.vector() %>% unlist() %>% as.vector() %>% mean() * 100
    
  } else{
    inlike[x] <- deduplicated_table %>% filter(Analysis.Region == "VI_IN_LIKE" & PanCK.Positive.Classification == 1) %>% dplyr::select(paste0(the_markers[x],".Positive.Classification")) %>% as.vector() %>% unlist() %>% as.vector() %>% mean() * 100
    dest[x] <- deduplicated_table %>% filter(Analysis.Region == "VI_DESTRUCTIVE" & PanCK.Positive.Classification == 1) %>% dplyr::select(paste0(the_markers[x],".Positive.Classification")) %>% as.vector() %>% unlist() %>% as.vector() %>% mean() * 100
    pdac[x] <- deduplicated_table %>% filter(Analysis.Region == "PDAC_enriched" & PanCK.Positive.Classification == 1) %>% dplyr::select(paste0(the_markers[x],".Positive.Classification")) %>% as.vector() %>% unlist() %>% as.vector() %>% mean() * 100
  }
}

# getting table of each percent positivity for each marker by tissue type (analysis region)
output_table <- as.data.frame(cbind(inlike, dest, pdac))

# making the row names the marker names
rownames(output_table) <- the_markers

# adding column for marker names
output_table$gene <- rownames(output_table)

# pivot long for plotting
pivoted <- as.data.frame(pivot_longer(output_table, cols = c(1,2,3)))

# for the odds ratios, pdac will be the dominator, but when comparing in like to dest, dest will be the denominator
inlike <- deduplicated_table %>% filter(Analysis.Region == "VI_IN_LIKE") %>% dplyr::select(paste0(the_markers[2],".Positive.Classification")) %>% as.vector() %>% unlist() %>% as.vector()
dest <- deduplicated_table %>% filter(Analysis.Region == "VI_DESTRUCTIVE") %>% dplyr::select(paste0(the_markers[2],".Positive.Classification")) %>% as.vector() %>% unlist() %>% as.vector()

# Getting outcome and
inlike_vs_dest <- table(
  Population = c(rep("inlike", length(inlike)), rep("dest", length(dest))),
  Outcome = c(inlike, dest)
)
theodds <- oddsratio(inlike_vs_dest, method = "fisher")
# fisher p value
theodds$p.value[2,2]
# odds ratio
theodds$measure[2,1]

# doing odds ratios
inlike_vs_dest_list <- list()
inlike_vs_pdac_list <- list()
dest_vs_pdac_list <- list()
for(x in 1:length(the_markers)){
  if(the_markers[x] == "CLDN18.2")
  {
    inlike <- deduplicated_table %>% filter(Analysis.Region == "VI_IN_LIKE") %>% dplyr::select(paste0(the_markers[x],".Positive.Cytoplasm.Classification")) %>% as.vector() %>% unlist() %>% as.vector() 
    dest <- deduplicated_table %>% filter(Analysis.Region == "VI_DESTRUCTIVE") %>% dplyr::select(paste0(the_markers[x],".Positive.Cytoplasm.Classification"))  %>% as.vector() %>% unlist() %>% as.vector() 
    pdac <- deduplicated_table %>% filter(Analysis.Region == "PDAC_enriched") %>% dplyr::select(paste0(the_markers[x],".Positive.Cytoplasm.Classification")) %>% as.vector() %>% unlist() %>% as.vector() 
  } else{
    inlike <- deduplicated_table %>% filter(Analysis.Region == "VI_IN_LIKE") %>% dplyr::select(paste0(the_markers[x],".Positive.Classification")) %>% as.vector() %>% unlist() %>% as.vector() 
    dest <- deduplicated_table %>% filter(Analysis.Region == "VI_DESTRUCTIVE") %>% dplyr::select(paste0(the_markers[x],".Positive.Classification")) %>% as.vector() %>% unlist() %>% as.vector() 
    pdac <- deduplicated_table %>% filter(Analysis.Region == "PDAC_enriched") %>% dplyr::select(paste0(the_markers[x],".Positive.Classification")) %>% as.vector() %>% unlist() %>% as.vector() 
  }
  # Create a contingency table
  inlike_vs_dest <- table(
    Population = c(rep("2_inlike", length(inlike)), rep("1_dest", length(dest))),
    Outcome    = c(inlike, dest)
  )
  inlike_vs_dest_list[[x]] <- oddsratio(inlike_vs_dest, method = "fisher")
  
  names(inlike_vs_dest_list)[x] <- the_markers[x]
  
  # Create a contingency table
  inlike_vs_pdac <- table(
    Population = c(rep("2_inlike", length(inlike)), rep("1_pdac", length(pdac))),
    Outcome    = c(inlike, pdac)
  )
  inlike_vs_pdac_list[[x]] <- oddsratio(inlike_vs_pdac, method = "fisher")
  names(inlike_vs_pdac_list)[x] <- the_markers[x]
  
  # Create a contingency table
  dest_vs_pdac <- table(
    Population = c(rep("2_dest", length(dest)), rep("1_pdac", length(pdac))),
    Outcome    = c(dest, pdac)
  )
  dest_vs_pdac_list[[x]] <- oddsratio(dest_vs_pdac, method = "fisher")
  names(dest_vs_pdac_list)[x] <- the_markers[x]
  
}

oddslist <- list()
for(x in 1:length(the_markers)){
  comparisons <- c("inlike_vs_dest", "inlike_vs_pdac", "dest_vs_pdac")
  group_1 <- c("inlike", "inlike", "dest")
  group_2 <- c("dest", "pdac", "pdac")
  odds_ratios <- c(inlike_vs_dest_list[[x]]$measure[2,1], inlike_vs_pdac_list[[x]]$measure[2,1], dest_vs_pdac_list[[x]]$measure[2,1])
  pvalues <- c(inlike_vs_dest_list[[x]]$p.value[2,2], inlike_vs_pdac_list[[x]]$p.value[2,2], dest_vs_pdac_list[[x]]$p.value[2,2])
  thegene <- the_markers[x]
  oddslist[[x]]<- data.frame(comparisons, odds_ratios, pvalues, thegene, group_1, group_2)
  names(oddslist)[x] <- the_markers[x]
}

oddstable <- do.call(rbind, oddslist)
oddstable$p.signif <- sapply(oddstable$pvalues, p_to_asterisk)


# FIGURE 8B ---- Bar plots showing the percent positivity of all cells within each tissue type for each marker across all 19 patients ####

# TFF1 
tmp <- pivoted %>% filter(gene == the_markers[2] & name %in% c("inlike", "dest", "pdac"))
tmp$name <- factor(tmp$name, levels = c("inlike", "dest", "pdac"))
data_for_p_val_manual <- oddstable %>% filter(thegene == the_markers[2]) %>% dplyr::select("group_1", "group_2", "p.signif", "odds_ratios")
data_for_p_val_manual$label <- paste0("p = ", data_for_p_val_manual$p.signif, ", ", "OR = ", signif(data_for_p_val_manual$odds_ratios, digits = 2))
print(data_for_p_val_manual)
data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "label")
colnames(data_for_p_val_manual) <- c("group1", "group2", "label")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)
yposition <- max(tmp[,"value"]) * 1.05
stat.test <- add_pvalue(datman, y.position = yposition, step.increase = 0.12, label = "label", hide.ns = F, label.size = 4, bracket.nudge.y = 0.1, hjust = 0.5)
theplot <- ggplot() + geom_col(data = tmp, aes(x = name, y = value, fill = name)) +
  labs(fill = "", x = "", y = "Percentage of positive cells across all patients (%)", title = the_markers[2], color = "") +
  theme(axis.text.x = element_text(color = "black",
                                   size = 10, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  scale_fill_manual(name = "", labels = c("VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#36DAFF", "#9C27B0", "#DF2E2E")) +
  scale_x_discrete(labels= c("VI-IN-Like", "VI-Destructive", "PDAC")) +
  stat.test
pdf(paste0("outputs/plots_tables_objects/phenocycler/",the_markers[2],".pdf"), width = 4, height = 6)
print(theplot)
dev.off()

# CLDN18.2
tmp <- pivoted %>% filter(gene == the_markers[3] & name %in% c("inlike", "dest", "pdac"))
tmp$name <- factor(tmp$name, levels = c("inlike", "dest", "pdac"))
data_for_p_val_manual <- oddstable %>% filter(thegene == the_markers[3]) %>% dplyr::select("group_1", "group_2", "p.signif", "odds_ratios")
data_for_p_val_manual$label <- paste0("p = ", data_for_p_val_manual$p.signif, ", ", "OR = ", signif(data_for_p_val_manual$odds_ratios, digits = 2))
print(data_for_p_val_manual)
data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "label")
colnames(data_for_p_val_manual) <- c("group1", "group2", "label")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)
yposition <- max(tmp[,"value"]) * 1.05
stat.test <- add_pvalue(datman, y.position = yposition, step.increase = 0.12, label = "label", hide.ns = F, label.size = 4, bracket.nudge.y = 0.1, hjust = 0.5)
theplot <- ggplot() + geom_col(data = tmp, aes(x = name, y = value, fill = name)) +
  labs(fill = "", x = "", y = "Percentage of positive cells across all patients (%)", title = the_markers[3], color = "") +
  theme(axis.text.x = element_text(color = "black",
                                   size = 10, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  scale_fill_manual(name = "", labels = c("VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#36DAFF", "#9C27B0", "#DF2E2E")) +
  scale_x_discrete(labels= c("VI-IN-Like", "VI-Destructive", "PDAC")) +
  stat.test
pdf(paste0("outputs/plots_tables_objects/phenocycler/",the_markers[3],".pdf"), width = 4, height = 6)
print(theplot)
dev.off()

# GATA6
tmp <- pivoted %>% filter(gene == the_markers[4] & name %in% c("inlike", "dest", "pdac"))
tmp$name <- factor(tmp$name, levels = c("inlike", "dest", "pdac"))
data_for_p_val_manual <- oddstable %>% filter(thegene == the_markers[4]) %>% dplyr::select("group_1", "group_2", "p.signif", "odds_ratios")
data_for_p_val_manual$label <- paste0("p = ", data_for_p_val_manual$p.signif, ", ", "OR = ", signif(data_for_p_val_manual$odds_ratios, digits = 2))
print(data_for_p_val_manual)
data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "label")
colnames(data_for_p_val_manual) <- c("group1", "group2", "label")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)
yposition <- max(tmp[,"value"]) * 1.05
stat.test <- add_pvalue(datman, y.position = yposition, step.increase = 0.12, label = "label", hide.ns = F, label.size = 4, bracket.nudge.y = 0.1, hjust = 0.5)
theplot <- ggplot() + geom_col(data = tmp, aes(x = name, y = value, fill = name)) +
  labs(fill = "", x = "", y = "Percentage of positive cells across all patients (%)", title = the_markers[4], color = "") +
  theme(axis.text.x = element_text(color = "black",
                                   size = 10, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  scale_fill_manual(name = "", labels = c("VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#36DAFF", "#9C27B0", "#DF2E2E")) +
  scale_x_discrete(labels= c("VI-IN-Like", "VI-Destructive", "PDAC")) +
  stat.test
pdf(paste0("outputs/plots_tables_objects/phenocycler/",the_markers[4],".pdf"), width = 4, height = 6)
print(theplot)
dev.off()

# KRT5
tmp <- pivoted %>% filter(gene == the_markers[5] & name %in% c("inlike", "dest", "pdac"))
tmp$name <- factor(tmp$name, levels = c("inlike", "dest", "pdac"))
data_for_p_val_manual <- oddstable %>% filter(thegene == the_markers[5]) %>% dplyr::select("group_1", "group_2", "p.signif", "odds_ratios")
data_for_p_val_manual$label <- paste0("p = ", data_for_p_val_manual$p.signif, ", ", "OR = ", signif(data_for_p_val_manual$odds_ratios, digits = 2))
print(data_for_p_val_manual)
data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "label")
colnames(data_for_p_val_manual) <- c("group1", "group2", "label")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)
yposition <- max(tmp[,"value"]) * 1.05
stat.test <- add_pvalue(datman, y.position = yposition, step.increase = 0.12, label = "label", hide.ns = F, label.size = 4, bracket.nudge.y = 0.1, hjust = 0.45)
theplot <- ggplot() + geom_col(data = tmp, aes(x = name, y = value, fill = name)) +
  labs(fill = "", x = "", y = "Percentage of positive cells across all patients (%)", title = the_markers[5], color = "") +
  theme(axis.text.x = element_text(color = "black",
                                   size = 10, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  scale_fill_manual(name = "", labels = c("VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#36DAFF", "#9C27B0", "#DF2E2E")) +
  scale_x_discrete(labels= c("VI-IN-Like", "VI-Destructive", "PDAC")) +
  stat.test
pdf(paste0("outputs/plots_tables_objects/phenocycler/",the_markers[5],".pdf"), width = 4, height = 6)
print(theplot)
dev.off()

# KRT17
tmp <- pivoted %>% filter(gene == the_markers[6] & name %in% c("inlike", "dest", "pdac"))
tmp$name <- factor(tmp$name, levels = c("inlike", "dest", "pdac"))
data_for_p_val_manual <- oddstable %>% filter(thegene == the_markers[6]) %>% dplyr::select("group_1", "group_2", "p.signif", "odds_ratios")
data_for_p_val_manual$label <- paste0("p = ", data_for_p_val_manual$p.signif, ", ", "OR = ", signif(data_for_p_val_manual$odds_ratios, digits = 2))
print(data_for_p_val_manual)
data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "label")
colnames(data_for_p_val_manual) <- c("group1", "group2", "label")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)
yposition <- max(tmp[,"value"]) * 1.05
stat.test <- add_pvalue(datman, y.position = yposition, step.increase = 0.12, label = "label", hide.ns = F, label.size = 4, bracket.nudge.y = 0.1, hjust = 0.5)
theplot <- ggplot() + geom_col(data = tmp, aes(x = name, y = value, fill = name)) +
  labs(fill = "", x = "", y = "Percentage of positive cells across all patients (%)", title = the_markers[6], color = "") +
  theme(axis.text.x = element_text(color = "black",
                                   size = 10, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  scale_fill_manual(name = "", labels = c("VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#36DAFF", "#9C27B0", "#DF2E2E")) +
  scale_x_discrete(labels= c("VI-IN-Like", "VI-Destructive", "PDAC")) +
  stat.test
pdf(paste0("outputs/plots_tables_objects/phenocycler/",the_markers[6],".pdf"), width = 4, height = 6)
print(theplot)
dev.off()

# S100A2
tmp <- pivoted %>% filter(gene == the_markers[7] & name %in% c("inlike", "dest", "pdac"))
tmp$name <- factor(tmp$name, levels = c("inlike", "dest", "pdac"))
data_for_p_val_manual <- oddstable %>% filter(thegene == the_markers[7]) %>% dplyr::select("group_1", "group_2", "p.signif", "odds_ratios")
data_for_p_val_manual$label <- paste0("p = ", data_for_p_val_manual$p.signif, ", ", "OR = ", signif(data_for_p_val_manual$odds_ratios, digits = 2))
print(data_for_p_val_manual)
data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "label")
colnames(data_for_p_val_manual) <- c("group1", "group2", "label")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)
yposition <- max(tmp[,"value"]) * 1.05
stat.test <- add_pvalue(datman, y.position = yposition, step.increase = 0.12, label = "label", hide.ns = F, label.size = 4, bracket.nudge.y = 0.1, hjust = 0.45)
theplot <- ggplot() + geom_col(data = tmp, aes(x = name, y = value, fill = name)) +
  labs(fill = "", x = "", y = "Percentage of positive cells across all patients (%)", title = the_markers[7], color = "") +
  theme(axis.text.x = element_text(color = "black",
                                   size = 10, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  scale_fill_manual(name = "", labels = c("VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#36DAFF", "#9C27B0", "#DF2E2E")) +
  scale_x_discrete(labels= c("VI-IN-Like", "VI-Destructive", "PDAC")) +
  stat.test
pdf(paste0("outputs/plots_tables_objects/phenocycler/",the_markers[7],".pdf"), width = 4, height = 6)
print(theplot)
dev.off()

# CEACAM5
tmp <- pivoted %>% filter(gene == the_markers[8] & name %in% c("inlike", "dest", "pdac"))
tmp$name <- factor(tmp$name, levels = c("inlike", "dest", "pdac"))
data_for_p_val_manual <- oddstable %>% filter(thegene == the_markers[8]) %>% dplyr::select("group_1", "group_2", "p.signif", "odds_ratios")
data_for_p_val_manual$label <- paste0("p = ", data_for_p_val_manual$p.signif, ", ", "OR = ", signif(data_for_p_val_manual$odds_ratios, digits = 2))
print(data_for_p_val_manual)
data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "label")
colnames(data_for_p_val_manual) <- c("group1", "group2", "label")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)
yposition <- max(tmp[,"value"]) * 1.05
stat.test <- add_pvalue(datman, y.position = yposition, step.increase = 0.12, label = "label", hide.ns = F, label.size = 4, bracket.nudge.y = 0.1, hjust = 0.5)
theplot <- ggplot() + geom_col(data = tmp, aes(x = name, y = value, fill = name)) +
  labs(fill = "", x = "", y = "Percentage of positive cells across all patients (%)", title = the_markers[8], color = "") +
  theme(axis.text.x = element_text(color = "black",
                                   size = 10, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  scale_fill_manual(name = "", labels = c("VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#36DAFF", "#9C27B0", "#DF2E2E")) +
  scale_x_discrete(labels= c("VI-IN-Like", "VI-Destructive", "PDAC")) +
  stat.test
pdf(paste0("outputs/plots_tables_objects/phenocycler/",the_markers[8],".pdf"), width = 4, height = 6)
print(theplot)
dev.off()

# MUC13
tmp <- pivoted %>% filter(gene == the_markers[9] & name %in% c("inlike", "dest", "pdac"))
tmp$name <- factor(tmp$name, levels = c("inlike", "dest", "pdac"))
data_for_p_val_manual <- oddstable %>% filter(thegene == the_markers[9]) %>% dplyr::select("group_1", "group_2", "p.signif", "odds_ratios")
data_for_p_val_manual$label <- paste0("p = ", data_for_p_val_manual$p.signif, ", ", "OR = ", signif(data_for_p_val_manual$odds_ratios, digits = 2))
print(data_for_p_val_manual)
data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "label")
colnames(data_for_p_val_manual) <- c("group1", "group2", "label")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)
yposition <- max(tmp[,"value"]) * 1.05
stat.test <- add_pvalue(datman, y.position = yposition, step.increase = 0.25, label = "label", hide.ns = F, label.size = 4, bracket.nudge.y = 0.1, hjust = 0.5)
theplot <- ggplot() + geom_col(data = tmp, aes(x = name, y = value, fill = name)) +
  labs(fill = "", x = "", y = "Percentage of positive cells across all patients (%)", title = the_markers[9], color = "") +
  theme(axis.text.x = element_text(color = "black",
                                   size = 10, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  scale_fill_manual(name = "", labels = c("VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#36DAFF", "#9C27B0", "#DF2E2E")) +
  scale_x_discrete(labels= c("VI-IN-Like", "VI-Destructive", "PDAC")) +
  stat.test
pdf(paste0("outputs/plots_tables_objects/phenocycler/",the_markers[9],".pdf"), width = 4, height = 6)
print(theplot)
dev.off()

# LAMC2
tmp <- pivoted %>% filter(gene == the_markers[10] & name %in% c("inlike", "dest", "pdac"))
tmp$name <- factor(tmp$name, levels = c("inlike", "dest", "pdac"))
data_for_p_val_manual <- oddstable %>% filter(thegene == the_markers[10]) %>% dplyr::select("group_1", "group_2", "p.signif", "odds_ratios")
data_for_p_val_manual$label <- paste0("p = ", data_for_p_val_manual$p.signif, ", ", "OR = ", signif(data_for_p_val_manual$odds_ratios, digits = 2))
data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "label")
colnames(data_for_p_val_manual) <- c("group1", "group2", "label")
data_for_p_val_manual <- data_for_p_val_manual[c(1,2,3),]
datman <- tibble(data_for_p_val_manual)
yposition <- max(tmp[,"value"]) * 1.05
stat.test <- add_pvalue(datman, y.position = yposition, step.increase = 0.3, label = "label", hide.ns = F, label.size = 4, bracket.nudge.y = 0.1, hjust = 0.5)
theplot <- ggplot() + geom_col(data = tmp, aes(x = name, y = value, fill = name)) +
  labs(fill = "", x = "", y = "Percentage of positive cells across all patients (%)", title = the_markers[10], color = "") +
  theme(axis.text.x = element_text(color = "black",
                                   size = 10, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  scale_fill_manual(name = "", labels = c("VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#36DAFF", "#9C27B0", "#DF2E2E")) +
  scale_x_discrete(labels= c("VI-IN-Like", "VI-Destructive", "PDAC")) +
  stat.test
pdf(paste0("outputs/plots_tables_objects/phenocycler/",the_markers[10],".pdf"), width = 4, height = 6)
print(theplot)
dev.off()

#####

# All cells classical, basal, double positive
deduplicated_table$classical <- 0
deduplicated_table$basal <- 0
deduplicated_table$moffitt <- 0

# If cell is positive for any of the classical markers, classical column == 1
deduplicated_table$classical[which(deduplicated_table$GATA6.Positive.Classification == 1)] <- 1
deduplicated_table$classical[which(deduplicated_table$CLDN18.2.Positive.Cytoplasm.Classification == 1)] <- 1
deduplicated_table$classical[which(deduplicated_table$TFF1.Positive.Classification == 1)] <- 1

# If cell is positive for any of the basal-like markers, basal column == 1
deduplicated_table$basal[which(deduplicated_table$KRT17.Positive.Classification == 1)] <- 1
deduplicated_table$basal[which(deduplicated_table$KRT5.Positive.Cytoplasm.Classification == 1)] <- 1
deduplicated_table$basal[which(deduplicated_table$S100A2.Positive.Classification == 1)] <- 1

# Initiating column for moffitt classification
deduplicated_table$moffitt <- 0

# Making Moffitt classifications based on classical/basal-like positivity
deduplicated_table$moffitt[which(deduplicated_table$classical == 1 & deduplicated_table$basal == 0)] <- "classical"
deduplicated_table$moffitt[which(deduplicated_table$classical == 0 & deduplicated_table$basal == 1)] <- "basal_like"
deduplicated_table$moffitt[which(deduplicated_table$classical == 1 & deduplicated_table$basal == 1)] <- "double_positive"
deduplicated_table$moffitt[which(deduplicated_table$classical == 0 & deduplicated_table$basal == 0)] <- "double_negative"

# Counting number of cells by their tissue type (VI subtype, PDAC) and their Moffitt classification
in_like_classical <- deduplicated_table %>% filter(Analysis.Region == "VI_IN_LIKE" & moffitt == "classical") %>% nrow()
in_like_basal <- deduplicated_table %>% filter(Analysis.Region == "VI_IN_LIKE" & moffitt == "basal_like") %>% nrow()
in_like_dp <- deduplicated_table %>% filter(Analysis.Region == "VI_IN_LIKE" & moffitt == "double_positive") %>% nrow()
in_like_dn <- deduplicated_table %>% filter(Analysis.Region == "VI_IN_LIKE" & moffitt == "double_negative") %>% nrow()
dest_classical <- deduplicated_table %>% filter(Analysis.Region == "VI_DESTRUCTIVE" & moffitt == "classical") %>% nrow()
dest_basal <- deduplicated_table %>% filter(Analysis.Region == "VI_DESTRUCTIVE" & moffitt == "basal_like") %>% nrow()
dest_dp <- deduplicated_table %>% filter(Analysis.Region == "VI_DESTRUCTIVE" & moffitt == "double_positive") %>% nrow()
dest_dn <- deduplicated_table %>% filter(Analysis.Region == "VI_DESTRUCTIVE" & moffitt == "double_negative") %>% nrow()
pdac_classical <- deduplicated_table %>% filter(Analysis.Region == "PDAC_enriched" & moffitt == "classical") %>% nrow()
pdac_basal <- deduplicated_table %>% filter(Analysis.Region == "PDAC_enriched" & moffitt == "basal_like") %>% nrow()
pdac_dp <- deduplicated_table %>% filter(Analysis.Region == "PDAC_enriched" & moffitt == "double_positive") %>% nrow()
pdac_dn <- deduplicated_table %>% filter(Analysis.Region == "PDAC_enriched" & moffitt == "double_negative") %>% nrow()

# Getting the percentage of cells by moffitt classification for each tissue type
tmpin <- c(in_like_classical, in_like_basal, in_like_dp, in_like_dn)
tmpin <- tmpin/sum(tmpin) * 100
tmpdest <- c(dest_classical, dest_basal, dest_dp, dest_dn)
tmpdest <- tmpdest/sum(tmpdest) * 100
tmppdac <- c(pdac_classical, pdac_basal, pdac_dp, pdac_dn)
tmppdac <- tmppdac/sum(tmppdac) * 100

# Making data frame of the percent of each Moffitt label by tissue type
classical_basal <- data.frame("VI_IN_LIKE" = tmpin,
                              "VI_DESTRUCTIVE" = tmpdest,
                              "PDAC_enriched" = tmppdac)
classical_basal$moffitt <- c("Classical", "Basal-like", "Co-expressor", "Double negative")
classical_basal$moffitt <- factor(classical_basal$moffitt, levels = c("Classical", "Basal-like", "Co-expressor", "Double negative"))

# Pivot longer for plotting
classical_basal <- as.data.frame(pivot_longer(classical_basal, cols = c(1,2,3)))
colnames(classical_basal) <- c("Moffitt", "Tissue_type", "Percent")
classical_basal$Tissue_type <- factor(classical_basal$Tissue_type, levels = c("VI_IN_LIKE", "VI_DESTRUCTIVE", "PDAC_enriched"))

# FIGURE 8C Percent composition by Moffitt classification of each tissue type across all 19 samples ####
pdf("outputs/plots_tables_objects/phenocycler/Moffitt_all.pdf", width = 5, height = 6)
ggplot(classical_basal) + geom_col(aes(x = Tissue_type, y = Percent, fill = Moffitt), position = "stack")+
  scale_fill_manual(values = c("#8fff5a", "#ffb75a", "#6c6e00", "grey")) +
  labs(fill = "", x = "", y = "Percent positive", title = "Moffitt Classification", color = "") +
  theme(axis.text.x = element_text(color = "black",
                                   size = 8, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18))
dev.off()
#####

# Making figure 8C except by individual patient
thetmplist <- list()
# Looping through each patient
for(x in 1:length(unique(deduplicated_table$Image.Location))){
  tmp_total <- deduplicated_table %>% filter(Image.Location == unique(deduplicated_table$Image.Location)[x])
  # Counting number of cells by their tissue type (VI subtype, PDAC) and their Moffitt classification
  tmpic <- tmp_total %>% filter(Analysis.Region == "VI_IN_LIKE" & moffitt == "classical") %>% nrow()
  tmpib <- tmp_total %>% filter(Analysis.Region == "VI_IN_LIKE" & moffitt == "basal_like") %>% nrow()
  tmpidp <- tmp_total %>% filter(Analysis.Region == "VI_IN_LIKE" & moffitt == "double_positive") %>% nrow()
  tmpidn <- tmp_total %>% filter(Analysis.Region == "VI_IN_LIKE" & moffitt == "double_negative") %>% nrow()
  
  tmpdc <- tmp_total %>% filter(Analysis.Region == "VI_DESTRUCTIVE" & moffitt == "classical") %>% nrow()
  tmpdb <- tmp_total %>% filter(Analysis.Region == "VI_DESTRUCTIVE" & moffitt == "basal_like") %>% nrow()
  tmpddp <- tmp_total %>% filter(Analysis.Region == "VI_DESTRUCTIVE" & moffitt == "double_positive") %>% nrow()
  tmpddn <- tmp_total %>% filter(Analysis.Region == "VI_DESTRUCTIVE" & moffitt == "double_negative") %>% nrow()
  
  tmppc <- tmp_total %>% filter(Analysis.Region == "PDAC_enriched" & moffitt == "classical") %>% nrow()
  tmppb <- tmp_total %>% filter(Analysis.Region == "PDAC_enriched" & moffitt == "basal_like") %>% nrow()
  tmppdp <- tmp_total %>% filter(Analysis.Region == "PDAC_enriched" & moffitt == "double_positive") %>% nrow()
  tmppdn <- tmp_total %>% filter(Analysis.Region == "PDAC_enriched" & moffitt == "double_negative") %>% nrow()
  
  # Getting the percentage of cells by moffitt classification for each tissue type
  the_number <- c(tmpic, tmpib, tmpidp, tmpidn, tmpdc, tmpdb, tmpddp, tmpddn, tmppc, tmppb, tmppdp, tmppdn)
  the_percent <- c(the_number[1:4]/sum(the_number[1:4]) * 100, the_number[5:8]/sum(the_number[5:8]) * 100, the_number[9:12]/sum(the_number[9:12]) * 100)
  the_moffitt <- rep(c("classical", "basal_like", "double_positive", "double_negative"), 3)
  the_tissue_type <- c(rep("VI_IN_LIKE", 4), rep("VI_DESTRUCTIVE", 4), rep("PDAC_enriched", 4))
  the_slide <- rep(unique(deduplicated_table$Image.Location)[x], 12)
  divided_by_pdac <- c(the_percent[1]/the_percent[9], the_percent[2]/the_percent[10], the_percent[3]/the_percent[11], the_percent[4]/the_percent[12],
                       the_percent[5]/the_percent[9], the_percent[6]/the_percent[10], the_percent[7]/the_percent[11], the_percent[8]/the_percent[12],
                       the_percent[9]/the_percent[9], the_percent[10]/the_percent[10], the_percent[11]/the_percent[11], the_percent[12]/the_percent[12])
  # Just counting the number of patients that have more/less/equal of each Moffitt subtype than PDAC
  compared_to_pdac <- c()
  compared_to_pdac[which(divided_by_pdac > 1)] <- "more"
  compared_to_pdac[which(divided_by_pdac < 1)] <- "less"
  compared_to_pdac[which(divided_by_pdac == 1)] <- "equal"
  thetmpdf <- data.frame(the_number, the_percent, the_moffitt, the_tissue_type, the_slide, divided_by_pdac, compared_to_pdac)
  thetmplist[[x]] <- thetmpdf
}

# Getting all of the individual patient data into a single data frame
cbound <- do.call(cbind, thetmplist)

# in_like.. shows that 4 patients have IN-LIKE which are less classical than PDAC, while 10 patients have IN-LIKE that are more classical
# also shows that 13 patients have IN-LIKE that are less basal-like than PDAC, 1 which is more basal like
which(cbound[1,] == "less") %>% length()
which(cbound[1,] == "more") %>% length()
which(cbound[2,] == "less") %>% length()
which(cbound[2,] == "more") %>% length()
which(cbound[3,] == "less") %>% length()
which(cbound[3,] == "more") %>% length()
which(cbound[4,] == "less") %>% length()
which(cbound[4,] == "more") %>% length()
# destructive.. shows that 8/14 patients have destructive that are less classical than PDAC, 6/14 are more classical. 8/14 are less basal, 6/10 are more basal
which(cbound[5,] == "less") %>% length()
which(cbound[5,] == "more") %>% length()
which(cbound[6,] == "less") %>% length()
which(cbound[6,] == "more") %>% length()
which(cbound[7,] == "less") %>% length()
which(cbound[7,] == "more") %>% length()
which(cbound[8,] == "less") %>% length()
which(cbound[8,] == "more") %>% length()

# making the classical basal intra-patient comparison figures
individual_classical_basal <- do.call(rbind, thetmplist)
individual_classical_basal$the_moffitt <- factor(individual_classical_basal$the_moffitt, levels = c("classical", "basal_like", "double_positive", "double_negative"))
individual_classical_basal$the_tissue_type <- factor(individual_classical_basal$the_tissue_type, levels = c("VI_IN_LIKE", "VI_DESTRUCTIVE", "PDAC_enriched"))

# FIGURE S12D ---- Percentage of cells in each Moffitt subtype by tissue type, measured in each patient ####
pdf("outputs/plots_tables_objects/phenocycler/individual_moffitt.pdf", height = 6, width = 6)
ggplot(individual_classical_basal) + geom_col(aes(x = the_tissue_type, y = the_percent, fill = the_moffitt), position = "stack")+
  facet_wrap(~the_slide, ncol =1, strip.position = "right") +
  theme(strip.background = element_blank(), strip.placement = "outside", strip.text.y = element_text(angle= 0)) +
  scale_fill_manual(c("Classical", "Basal-like", "Co-expressor", "Double negative"), values = c("#8fff5a", "#ffb75a", "#6c6e00", "grey")) +
  labs(fill = "", x = "", y = "Percent positive", title = "Moffitt Classification", color = "") +
  theme(axis.text.x = element_text(color = "black",
                                   size = 8, angle = 45, hjust = 1), axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18))
dev.off()
#####

# All foci classical, basal, double positive. double positive foci = at least 25% of the foci is classical (or basal) or double positive AKA a single phenotype is not found in more than 50% of the cells in the foci

# Initializing classical and basal-like columns
annotated_total
annotated_total$classical <- 0
annotated_total$basal <- 0

# If cell is positive for any classical marker, classical column == 1 
annotated_total$classical[which(annotated_total$GATA6.Positive.Classification == 1)] <- 1
annotated_total$classical[which(annotated_total$CLDN18.2.Positive.Cytoplasm.Classification == 1)] <- 1
annotated_total$classical[which(annotated_total$TFF1.Positive.Classification == 1)] <- 1

# If cell is positive for any basal-like marker, basal column == 1 
annotated_total$basal[which(annotated_total$KRT17.Positive.Classification == 1)] <- 1
annotated_total$basal[which(annotated_total$KRT5.Positive.Classification == 1)] <- 1
annotated_total$basal[which(annotated_total$S100A2.Positive.Classification == 1)] <- 1

# Initializing moffitt classification column
annotated_total$moffitt <- 0

# Classical and basal-like groups
annotated_total$moffitt[which(annotated_total$classical == 1 & annotated_total$basal == 0)] <- "classical"
annotated_total$moffitt[which(annotated_total$classical == 0 & annotated_total$basal == 1)] <- "basal_like"
annotated_total$moffitt[which(annotated_total$classical == 1 & annotated_total$basal == 1)] <- "double_positive"
annotated_total$moffitt[which(annotated_total$classical == 0 & annotated_total$basal == 0)] <- "double_negative"
annotated_total$moffitt_foci <- 0
annotated_total %>% filter(ROI_by_slide == unique(ROI_by_slide)[1])

# Looping through each patient and creating the Moffitt classifications for each foci based on their representation of classical/basal markers
for(x in 1:length(unique(annotated_total$ROI_by_slide))){
  tmp <- annotated_total %>% filter(ROI_by_slide == unique(ROI_by_slide)[x])
  tmpc <- mean(tmp$classical)
  tmpb <- mean(tmp$basal)
  if(tmpb >= 0.25 & tmpc >= 0.25){
    annotated_total$moffitt_foci[which(annotated_total$ROI_by_slide %in% tmp$ROI_by_slide)] <- "double_positive"
  }
  else{
    annotated_total$moffitt_foci[which(annotated_total$ROI_by_slide %in% tmp$ROI_by_slide)] <- "double_negative"
  }
  if(tmpc > 0.50 & tmpb <= 0.25){
    annotated_total$moffitt_foci[which(annotated_total$ROI_by_slide %in% tmp$ROI_by_slide)] <- "classical"
  }
  if(tmpb > 0.50 & tmpc <= 0.25){
    annotated_total$moffitt_foci[which(annotated_total$ROI_by_slide %in% tmp$ROI_by_slide)] <- "basal_like"
  }
  print(paste0("Loop ",x," completed"))
}

# Getting vector of names for the individual ROIs (foci)
unique_roi <- annotated_total$ROI_by_slide[!duplicated(annotated_total$ROI_by_slide)]

# Getting the corresponding tissue type for each of the above ROIs (foci)
tissue_by_roi <- annotated_total$Analysis.Region[!duplicated(annotated_total$ROI_by_slide)]

# Getting the corresponding moffitt classification for each of the ROIs (foci)
moffitt_by_roi <- annotated_total$moffitt_foci[!duplicated(annotated_total$ROI_by_slide)]

# Making data frame with all ROIs, their tissue type and their Moffitt group
thedf <- data.frame(unique_roi, tissue_by_roi, moffitt_by_roi)
colnames(thedf) <- c("ROI_by_slide", "Analysis.Region", "moffitt_foci")

# getting the number of foci in each moffitt subtype by tissue type
infoci <- thedf %>% filter(Analysis.Region == "VI_IN_LIKE") %>% dplyr::select("moffitt_foci") %>% as.vector() %>% unlist() %>% as.vector() %>% table()
destfoci <- thedf %>% filter(Analysis.Region == "VI_DESTRUCTIVE") %>% dplyr::select("moffitt_foci") %>% as.vector() %>% unlist() %>% as.vector() %>% table()
pdacfoci <- thedf %>% filter(Analysis.Region == "PDAC_enriched") %>% dplyr::select("moffitt_foci") %>% as.vector() %>% unlist() %>% as.vector() %>% table()

# putting zero here because there are zero basal-like IN-Like foci
infoci[4] <- 0
names(infoci)[4] <- "basal_like"

# switching the column order so that it is consistent with the other three
infoci <- c(infoci[4], infoci[1:3])

# Getting percent of foci in each Moffitt subtype by each tissue type 
percentinfoci <- infoci/sum(infoci) * 100
percentdestfoci <- destfoci/sum(destfoci) * 100
percentpdacfoci <- pdacfoci/sum(pdacfoci) * 100

# Making a data frame with the number and percent of each Moffitt subtype by foci of each tissue type
classical_basal_foci <- data.frame(moffitt = rep(names(destfoci),3),  
                                   tissue_type = c(rep("VI_IN_LIKE", 4), rep("VI_DESTRUCTIVE", 4), rep("PDAC", 4)),
                                   tissue_count = c(infoci, destfoci, pdacfoci),
                                   tissue_percent = c(percentinfoci, percentdestfoci, percentpdacfoci))

# factoring by tissue type and moffitt classification for plotting 
classical_basal_foci$tissue_type <- factor(classical_basal_foci$tissue_type, levels =c("VI_IN_LIKE", "VI_DESTRUCTIVE", "PDAC"))
classical_basal_foci$moffitt <- factor(classical_basal_foci$moffitt, levels =c("classical", "basal_like", "double_positive", "double_negative"))

# FIGURE S12C ---- Percent composition of each tissue type across all samples by Moffitt classification of individual foci ####
pdf("outputs/plots_tables_objects/phenocycler/Moffitt_all_by_foci.pdf", width = 5, height = 6)
ggplot(classical_basal_foci) + geom_col(aes(x = tissue_type, y = tissue_percent, fill = moffitt), position = "stack")+
  scale_fill_manual(values = c("#8fff5a", "#ffb75a", "#6c6e00", "grey")) +
  labs(fill = "", x = "", y = "Percent positive", title = "Moffitt Classification", color = "") +
  theme(axis.text.x = element_text(color = "black",
                                   size = 8, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18))
dev.off()
#####

# Each marker's percent positivity within each patient, made for between-patient comparisons. This code compares the per-patient percent positivity, which controls for cell number in each patient
output_by_slide <- list()
# Looping through each patient, getting the percent positivity for each marker, for each tissue type, for each patient
for(x in 1:length(unique(deduplicated_table$Image.Location))){
  theimage <- unique(deduplicated_table$Image.Location)[x]
  theindex <- which(deduplicated_table$Image.Location == theimage)
  tmp_table <- deduplicated_table[theindex,]
  inlike <- c()
  dest <- c()
  pdac <- c()
  for(y in 1:length(the_markers)){
    if(the_markers[y] == "CLDN18.2")
    {
      inlike[y] <- tmp_table %>% filter(Analysis.Region == "VI_IN_LIKE" & PanCK.Positive.Classification == 1) %>% dplyr::select(paste0(the_markers[y],".Positive.Cytoplasm.Classification")) %>% as.vector() %>% unlist() %>% as.vector() %>% mean() * 100
      dest[y] <- tmp_table %>% filter(Analysis.Region == "VI_DESTRUCTIVE" & PanCK.Positive.Classification == 1) %>% dplyr::select(paste0(the_markers[y],".Positive.Cytoplasm.Classification")) %>% as.vector() %>% unlist() %>% as.vector() %>% mean() * 100
      pdac[y] <- tmp_table %>% filter(Analysis.Region == "PDAC_enriched" & PanCK.Positive.Classification == 1) %>% dplyr::select(paste0(the_markers[y],".Positive.Cytoplasm.Classification")) %>% as.vector() %>% unlist() %>% as.vector() %>% mean() * 100
      
    } else{
      inlike[y] <- tmp_table %>% filter(Analysis.Region == "VI_IN_LIKE" & PanCK.Positive.Classification == 1) %>% dplyr::select(paste0(the_markers[y],".Positive.Classification")) %>% as.vector() %>% unlist() %>% as.vector() %>% mean() * 100
      dest[y] <- tmp_table %>% filter(Analysis.Region == "VI_DESTRUCTIVE" & PanCK.Positive.Classification == 1) %>% dplyr::select(paste0(the_markers[y],".Positive.Classification")) %>% as.vector() %>% unlist() %>% as.vector() %>% mean() * 100
      pdac[y] <- tmp_table %>% filter(Analysis.Region == "PDAC_enriched" & PanCK.Positive.Classification == 1) %>% dplyr::select(paste0(the_markers[y],".Positive.Classification")) %>% as.vector() %>% unlist() %>% as.vector() %>% mean() * 100
    }
  }
  output_table <- as.data.frame(cbind(inlike, dest, pdac))
  rownames(output_table) <- the_markers
  output_table$slide <- unique(deduplicated_table$Image.Location)[x]
  output_by_slide[[x]] <- output_table
}

# Organizing the output into tables that show the percent positivity for each marker for each patient, one table for each tissue type
tmp <- as.data.frame(do.call(cbind, output_by_slide))
inlike_table <- tmp[,grep("inlike",colnames(tmp))]
dest_table <- tmp[,grep("dest",colnames(tmp))]
pdac_table <- tmp[,grep("pdac",colnames(tmp))]

# removing NA (patients that don't contain the tissue type measured in the table)
inlike_table <- inlike_table[,colSums(is.na(inlike_table))<nrow(inlike_table)]
dest_table <- dest_table[,colSums(is.na(dest_table))<nrow(dest_table)]
pdac_table <- pdac_table[,colSums(is.na(pdac_table))<nrow(pdac_table)]

# pivot longer 
inlike_table$genes <- rownames(inlike_table)
inlike_table <- pivot_longer(inlike_table, cols = c(1:14)) %>% as.data.frame()
inlike_table$name <- "VI_IN_LIKE"

# pivot longer
dest_table$genes <- rownames(dest_table)
thecol <- ncol(dest_table)-1
dest_table <- pivot_longer(dest_table, cols = c(1:thecol)) %>% as.data.frame()
dest_table$name <- "VI_DESTRUCTIVE"

# pivot longer
pdac_table$genes <- rownames(pdac_table)
thecol <- ncol(pdac_table)-1
pdac_table <- pivot_longer(pdac_table, cols = c(1:thecol)) %>% as.data.frame()
pdac_table$name <- "PDAC"

# Make table that has column for genes, tissue type and percent positivity for each of the genes
the_averaged_table <- rbind(inlike_table, dest_table, pdac_table)
the_averaged_table$name <- factor(the_averaged_table$name, levels = c("VI_IN_LIKE", "VI_DESTRUCTIVE", "PDAC"))
the_averaged_table$genes <- factor(the_averaged_table$genes, levels = c("TFF1", "GATA6", "CLDN18.2", "KRT17", "KRT5", "S100A2", "CEACAM5", "MUC13", "LAMC2"))

# FIGURE S12B ---- Box-plot showing the per-patient proportion of cell positivity for each marker by tissue type ####
pdf("outputs/plots_tables_objects/phenocycler/Average_of_individual_slides.pdf")
ggplot(data = subset(the_averaged_table, genes %in% c("TFF1", "GATA6", "CLDN18.2", "KRT17", "KRT5", "S100A2", "CEACAM5", "MUC13", "LAMC2"))) + geom_boxplot(aes(x = genes, y = value, fill = name)) +
  scale_fill_manual(name = "", labels = c("VI-IN-Like", "VI-Destructive", "PDAC"), values = c("#36DAFF", "#9C27B0", "#DF2E2E")) +
  labs(fill = "", x = "", y = "Percent positive", title = "Averages of individual samples", color = "") +
  theme(axis.text.x = element_text(color = "black",
                                   size = 12, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18))
dev.off()
#####

# This gives me all of the cells which did not belong to a tile with at least 10 cells
the_missing_cells <- deduplicated_table[which(deduplicated_table$thecellnames %notin% annotated_total$thecellnames),]
the_missing_cells$ROI <- NA
the_missing_cells$ROI_by_slide <- NA
the_missing_cells$moffitt_foci <- NA

the_missing_cells <- the_missing_cells[,colnames(annotated_total)]
fully_annotated_panck_positive <- rbind(annotated_total, the_missing_cells)

# fully annotated epithelial phenocycler data
write.csv(fully_annotated_panck_positive, "outputs/plots_tables_objects/phenocycler/processed_phenocycler_data.csv")

# ANNOTATING THE VI FOCI

# Here, we used a custom shiny app to manually annotate individual VI foci, since HALO was not able to import annotations for individual foci.
# Each tissue section was manually annotated for each individual VI foci/subtype
# Manually change value of x (1-19) and then open shiny app to view manual annotation
x <- 1
for_annotation <- deduplicated_table %>% filter(Image.Location == unique(deduplicated_table$Image.Location)[x] & Analysis.Region %in% c("PNI", "VI_IN_LIKE", "VI_DESTRUCTIVE", "VI_CONVENTIONAL", "Possible_PNI", "NORMAL_DUCT"))
for_annotation_pos <- for_annotation %>% filter(PanCK.Positive.Classification == 1)
for_annotation_pos$names <- rownames(for_annotation_pos)
for_annotation_pos_coords <- for_annotation_pos %>% dplyr::select("XMin", "XMax", "YMin", "YMax", "Analysis.Region")
for_annotation_pos_coords$xmean <- (for_annotation_pos_coords$XMin + for_annotation_pos_coords$XMax)/2
for_annotation_pos_coords$ymean <- (for_annotation_pos_coords$YMin + for_annotation_pos_coords$YMax)/2
for_annotation_pos_coords <- for_annotation_pos_coords %>% dplyr::select("xmean", "ymean", "Analysis.Region")
for_annotation_pos_coords$ymean <- for_annotation_pos_coords$ymean * -1

# Define the global variable in the global environment
global_named_points <- list()

# Define UI
ui <- fluidPage(
  titlePanel("Select Points by Dragging a Box, Name Them, and Save"),
  sidebarLayout(
    sidebarPanel(
      h4("Selected Points"),
      verbatimTextOutput("selected_points"),
      textInput("point_name", "Enter a name for the selected points:"),
      actionButton("save_button", "Save Named Points"),
      br(), br(),
      h4("Saved Named Points"),
      verbatimTextOutput("saved_points")
    ),
    mainPanel(
      plotOutput("scatter_plot", brush = brushOpts(id = "plot_brush"))
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Create a dataset
  set.seed(123)
  data <- for_annotation_pos_coords
  
  # Reactive value to store selected points
  selected_points <- reactiveVal(NULL)
  
  # Render the scatter plot
  output$scatter_plot <- renderPlot({
    ggplot(data, aes(x = xmean, y = ymean, col = Analysis.Region)) +
      geom_point(size = 3) +
      theme_minimal()
  })
  
  # Update selected points when brushing
  observeEvent(input$plot_brush, {
    brushed <- brushedPoints(data, input$plot_brush, xvar = "xmean", yvar = "ymean")
    selected_points(brushed)
  })
  
  # Display selected points
  output$selected_points <- renderPrint({
    selected_points()
  })
  
  # Save named points to the global variable
  observeEvent(input$save_button, {
    if (!is.null(selected_points())) {
      name <- input$point_name
      if (name == "") {
        showModal(modalDialog(
          title = "Error",
          "Please enter a name for the selected points.",
          easyClose = TRUE
        ))
      } else {
        # Attach the name to the selected points
        named_points <- selected_points()
        named_points$name <- name
        
        # Save to the global variable
        global_named_points[[name]] <<- named_points  # Use <<- to assign to the global environment
        
        # Clear the selected points and name input
        selected_points(NULL)
        updateTextInput(session, "point_name", value = "")
        
        showModal(modalDialog(
          title = "Success",
          paste("Named points have been saved as:", name),
          easyClose = TRUE
        ))
      }
    } else {
      showModal(modalDialog(
        title = "Error",
        "No points selected.",
        easyClose = TRUE
      ))
    }
  })
  
  # Display saved named points
  output$saved_points <- renderPrint({
    global_named_points
  })
}

# Run the app
shinyApp(ui = ui, server = server)
names(global_named_points)
#codex_01_groups <- global_named_points
#saveRDS(codex_01_groups, "outputs/plots_tables_objects/phenocycler/codex_01_groups.rds")

#codex_05_groups <- global_named_points
#saveRDS(codex_05_groups, "outputs/plots_tables_objects/phenocycler/codex_05_groups.rds")

#codex_09_groups <- global_named_points
#saveRDS(codex_09_groups, "outputs/plots_tables_objects/phenocycler/codex_09_groups.rds")

#codex_10_groups <- global_named_points
#saveRDS(codex_10_groups, "outputs/plots_tables_objects/phenocycler/codex_10_groups.rds")

#codex_12_groups <- global_named_points
#saveRDS(codex_12_groups, "outputs/plots_tables_objects/phenocycler/codex_12_groups.rds")

#codex_14_groups <- global_named_points
#saveRDS(codex_14_groups, "outputs/plots_tables_objects/phenocycler/codex_14_groups.rds")

#codex_16_groups <- global_named_points
#saveRDS(codex_16_groups, "outputs/plots_tables_objects/phenocycler/codex_16_groups.rds")

#codex_17_groups <- global_named_points
#saveRDS(codex_17_groups, "outputs/plots_tables_objects/phenocycler/codex_17_groups.rds")

#codex_19_groups <- global_named_points
#saveRDS(codex_19_groups, "outputs/plots_tables_objects/phenocycler/codex_19_groups.rds")

# codex_21_groups <- global_named_points
# codex_21$IN_Like_1 <- codex_21$IN_Like_1 %>% filter(Analysis.Region == "VI_IN_LIKE")
# codex_21$Conventional_2 <- codex_21$Conventional_2 %>% filter(Analysis.Region == "VI_CONVENTIONAL")
# saveRDS(codex_21, "outputs/plots_tables_objects/phenocycler/codex_21_groups.rds")

# codex_21_groups$IN_Like_1 <- codex_21_groups$IN_Like_1 %>% filter(Analysis.Region == "VI_IN_LIKE")
# codex_21_groups$Conventional_2 <- codex_21_groups$Conventional_2 %>% filter(Analysis.Region == "VI_CONVENTIONAL")
# saveRDS(codex_21_groups, "outputs/plots_tables_objects/phenocycler/codex_21_groups.rds")

#codex_23_groups <- global_named_points
#saveRDS(codex_23_groups, "outputs/plots_tables_objects/phenocycler/codex_23_groups.rds")

#codex_25_groups <- global_named_points
#saveRDS(codex_25_groups, "outputs/plots_tables_objects/phenocycler/codex_25_groups.rds")

#codex_26_groups <- global_named_points
#saveRDS(codex_26_groups, "outputs/plots_tables_objects/phenocycler/codex_26_groups.rds")

#codex_28_groups <- global_named_points
#saveRDS(codex_28_groups, "outputs/plots_tables_objects/phenocycler/codex_28_groups.rds")

#codex_29_groups <- global_named_points
#saveRDS(codex_29_groups, "outputs/plots_tables_objects/phenocycler/codex_29_groups.rds")

#codex_32_groups <- global_named_points
#saveRDS(codex_32_groups, "outputs/plots_tables_objects/phenocycler/codex_32_groups.rds")

#KH_002_groups <- global_named_points
#saveRDS(KH_002_groups, "outputs/plots_tables_objects/phenocycler/KH_002_groups.rds")

#KH_003_groups <- global_named_points
#saveRDS(KH_003_groups, "outputs/plots_tables_objects/phenocycler/KH_003_groups.rds")

#KH_005_groups <- global_named_points
#saveRDS(KH_005_groups, "outputs/plots_tables_objects/phenocycler/KH_005_groups.rds")