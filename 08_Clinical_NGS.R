library(tidyverse)
library(reshape2)
library(limma)
library(data.table)
library(ggrepel)
library(ggpubr)
library(ggplot2)
library(ComplexHeatmap)
library(gridExtra)
library(cowplot)
library(psych)
library(FSA)
library(sva)
library(Biobase)
library(patchwork)
library(DESeq2)
library(caroline)
library('DT')
library('EnhancedVolcano')
library(webshot)
library(VennDiagram)
library(rstatix)
library(projectR)
library(fgsea)
library(GSVA)
library(msigdbr)
sessionInfo()

source("scripts/00_Custom_Functions.R")

# importing the genetic mutations from the NGS assay
vec_gene_238 <- c("HRAS", "KRAS", "NF1", "PIK3CA", "TP53", "AXL", "BRD4", "CDH1", "DDX41", "FANCA", "FOXO1", "KDM4C", "MYO18A", "PAG1", "PLCG2", "RNF43", "TGFBR2", "WDR90")
vec_sample_238 <- rep("pt238", length(vec_gene_238))
vec_base_change_238 <- c("G>A", "C>A", "C>T", "C>A", "G>A", "A>C", "C>G", "T>C", "T>C", "C>G", "C>G", "A>C", "G>A", "G>A", "A>C", "G>A", "T>G", "T>G")
vec_aa_change_238 <- c("p.T13I", "p.G12V", "p.Q2218X", "p.Q60K", "p.R213X", "p.Q361P", "p.E1270D", "p.N1122S", "p.I180V", "p.C625S", "p.A99P", "p.N789T", "p.R1354W", "p.P46L", "p.Q387P", "p.R127W", "p.Y473D", "p.F758C")
vec_vaf_238 <- c(48.98, 10.67, 14.24, 44.53, 9.23, 46.22, 52.11, 47.1, 50.15, 43.33, 5.06, 40.95, 43.15, 45.97, 46.81, 9.86, 17.92, 52.06)
vec_tmr_238 <- rep(NA, length(vec_gene_238))
table_238 <- as.data.frame(cbind(vec_sample_238, vec_gene_238, vec_base_change_238, vec_aa_change_238, vec_vaf_238, vec_tmr_238))

vec_gene_2N <- c("KRAS", "TP53", "TSC2", "CUX1", "ELP2", "EML4", "H3F3A", "KMT2B", "SMAD4")
vec_sample_2N <- rep("2N", length(vec_gene_2N))
vec_base_change_2N <- c("C>G", "C>A", "G>A", "G>A", "G>C", "A>G", "G>A", "A>G", "A>G")
vec_aa_change_2N <- c("p.G12R", "p.V272L", "p.G440S", "p.R570Q", "p.S491T", "p.H508R", "p.A76T", "p.K376E", "p.Y513C")
vec_vaf_2N <- c(45.74, 35.84, 48.24, 69.82, 67.07, 25.92, 20.09, 12.41, 16.45)
vec_tmr_2N <- rep(NA, length(vec_gene_2N))
table_2N <- as.data.frame(cbind(vec_sample_2N, vec_gene_2N, vec_base_change_2N, vec_aa_change_2N, vec_vaf_2N, vec_tmr_2N))

vec_gene_3K <- c("KRAS", "TP53", "CREBBP", "ECT2L", "IRF1", "PIK3CB", "TGFBR2")
vec_sample_3K <- rep("3K", length(vec_gene_3K))
vec_base_change_3K <- c("C>A", "C>A", "C>A", "G>A", "A>G", "C>T", "G>A")
vec_aa_change_3K <- c("p.G12V", "p.G245V", "p.Q946H", "p.R711Q", "p.V175A", "p.V199I", "p.V412M")
vec_vaf_3K <- c(9.95, 11.71, 44.45, 45.98, 35.91, 48.58, 60.79)
vec_tmr_3K <- rep(NA, length(vec_gene_3K))
table_3K <- as.data.frame(cbind(vec_sample_3K, vec_gene_3K, vec_base_change_3K, vec_aa_change_3K, vec_vaf_3K, vec_tmr_3K))

vec_gene_499 <- c("ATM","KRAS","BCL2","BCR", "BRD4", "BRD4", "CHD1", "NOTCH1")
vec_sample_499 <- rep("pt499", length(vec_gene_499))
vec_base_change_499 <- c("T>C", "C>A", "G>C","C>T","CGCTGCT>C","G>T", "C>T", "C>T")
vec_aa_change_499 <- c("p.V410A", "p.G12V", "p.A32G","p.R1080C", "p.Q1288_Q1289del", "p.P1019H", "p.D133N", "p.E848K")
vec_vaf_499 <- c(44.86, 6.32, 51, 9.73, 35.68, 55.95, 45.36, 49.84)
vec_tmr_499 <- rep(NA, length(vec_gene_499))
table_499 <- as.data.frame(cbind(vec_sample_499, vec_gene_499, vec_base_change_499, vec_aa_change_499, vec_vaf_499, vec_tmr_499))

vec_gene_602 <- c("CDKN2A", "KRAS", "SMO", "STK11", "TP53", "ACTB", "ARID1B", "ARID1B", "FOXO1", "GNA13", "HSP90AA1", "MKI67", "NSD1", "NSD1", "SMAD4", "TET2")
vec_sample_602 <- rep("pt602", length(vec_gene_602))
vec_base_change_602 <- c("G>A", "C>T", "C>G", "G>A", "C>T", "G>A", "G>GAGGAGC", "G>C", "C>G", "C>G", "C>T", "G>A","C>G","G>C","C>T","C>T")
vec_aa_change_602 <- c("p.R80X", "p.G12D", "p.F252L", "p.D343N", "p.R273H", "UTR5", "p.G332_G333insAG", "p.G333A", "p.A99P", "p.V221L", "p.E608K", "p.P1964S", "p.L1317V", "p.E1332Q", "p.R135X", "p.S460F")
vec_vaf_602 <- c(38.29, 25.38, 52.24, 50.13, 33.03, 51.86, 18.21, 5.85, 5.28, 53.78, 28.62,45.48, 45.08, 45.83, 34.48, 46.3)
vec_tmr_602 <- rep(NA, length(vec_gene_602))
table_602 <- as.data.frame(cbind(vec_sample_602, vec_gene_602, vec_base_change_602, vec_aa_change_602, vec_vaf_602, vec_tmr_602))

first_vec_gene_2M <- c("KRAS", "FBXW7", "ASMTL", "ATM", "CD274", "EML4", "FHIT", "JAK3", "KAT6A", "KMT2B", "MAP2K4A", "POLE", "PTPRO", "SMAD4", "SNYE1", "SNYE1", "WDR90")
first_vec_sample_2M<- rep("2M_f", length(first_vec_gene_2M))
first_vec_base_change_2M <- c("T>G", "C>CAGG", "T>G", "A>T", "G>C", "C>G", "C>T", "G>A", "T>C", "G>A", "A>C", "AG>A", "A>G", "C>CG", "T>G", "C>T", "A>G")
first_vec_aa_change_2M <- c("p.Q61H","p.T15_G16insP", "p.E468A", "p.D1853V", "p.M36I", "p.S821C", "c.103+1G>A", "p.R840C", "p.Q1091R", "p.R2092Q", "p.S16R", "p.267Lfs*27", "p.H362R", "p.A532Gfs*44", "p.I7014L", "p.R4420Q", "p.S1421G")
first_vec_vaf_2M <- c(36.18, 30.15, 48.29, 49.15, 49.46, 56.96, 5.44, 48.82, 40.21, 49.38, 47.02, 46.64, 49.32, 20.68, 56.77, 55.82, 49.14)
first_vec_tmr_2M <- rep(1.76, length(first_vec_gene_2M))
table_2M_f <- as.data.frame(cbind(first_vec_sample_2M, first_vec_gene_2M, first_vec_base_change_2M, first_vec_aa_change_2M, first_vec_vaf_2M, first_vec_tmr_2M))

first_vec_gene_3K <- c("CREBBP", "ECT2L", "IRF1", "KRAS", "PIK3CB", "TGFBR2", "TP53")
first_vec_sample_3K<- rep("3K_f", length(first_vec_gene_3K))
first_vec_base_change_3K <- c("C>A", "G>A", "A>G", "C>A", "C>T", "G>A", "C>A")
first_vec_aa_change_3K <- c("p.Q946H", "p.R711Q", "p.V175A", "p.G12V", "p.V199I", "p.V412M", "p.G245V")
first_vec_vaf_3K <- c(45.95, 51.27, 43.02, 8.75, 44.31, 57.86, 8.04)
first_vec_tmr_3K <- rep(0.88, length(first_vec_gene_3K))
table_3K_f <- as.data.frame(cbind(first_vec_sample_3K, first_vec_gene_3K, first_vec_base_change_3K, first_vec_aa_change_3K, first_vec_vaf_3K, first_vec_tmr_3K))

first_vec_gene_4L <- c("BRIP1", "FGFR3", "KIT", "KRAS", "TP53", "U2AF1", "BARD1", "GRIN2A", "IRF8", "JAK2", "KMT2D", "TRAF3", "WDR90")
first_vec_sample_4L<- rep("4L_f", length(first_vec_gene_4L))
first_vec_base_change_4L <- c("G>T", "T>C", "A>G", "C>A", "T>C", "G>T", "T>C", "T>A", "A>C", "G>A", "T>TC", "C>T", "G>A")
first_vec_aa_change_4L <- c("p.S624*", "p.F384L", "p.Y568C", "p.G12V", "p.N239D", "p.S34Y", "p.E223G", "p.1043S", "p.T134P", "p.R1063H", "p.E1287fs", "p.R118W", "p.R578H")
first_vec_vaf_4L <- c(46.4, 48.82, 18.85, 15.45, 18.25, 13.91, 47.25, 49.46, 18.84, 57.03, 15.24, 47.98, 51.24)
first_vec_tmr_4L <- rep(4.40, length(first_vec_gene_4L))
table_4L_f <- as.data.frame(cbind(first_vec_sample_4L, first_vec_gene_4L, first_vec_base_change_4L, first_vec_aa_change_4L, first_vec_vaf_4L, first_vec_tmr_4L))

first_vec_gene_499 <- c("ATM", "KRAS", "BCL2", "BCR", "EXOSC6", "NOTCH1")
first_vec_sample_499 <- rep("pt499_f", length(first_vec_gene_499))
first_vec_base_change_499 <- c("T>C", "C>A", "G>C", "C>T", "G>A", "C>T")
first_vec_aa_change_499 <- c("p.V410A", "p.G12V", "p.A32G", "p.R1080C", "p.P2L", "p.E848K")
first_vec_vaf_499 <- c(45.18, 4.92, 47.38, 11.13, 47.97, 47.5)
first_vec_tmr_499 <- rep(2.64, length(first_vec_gene_499))
table_499_f <- as.data.frame(cbind(first_vec_sample_499, first_vec_gene_499, first_vec_base_change_499, first_vec_aa_change_499, first_vec_vaf_499, first_vec_tmr_499))

first_vec_gene_063 <- c("ROS1", "CDKN2A", "KRAS", "TP53", "BLM", "HDAC7", "JAK2", "NOD1", "PDCD1LG2", "RECQL4")
first_vec_sample_063 <- rep("pt063_f", length(first_vec_gene_063))
first_vec_base_change_063 <- c("T>C", "G>A", "C>T", "G>A", "T>C", "G>A", "G>C", "C>T", "G>C", "C>T")
first_vec_aa_change_063 <- c("p.M633V", "p.H83Y", "p.G12D", "p.R196*", "p.V4A", "p.P300S", "p.E846D", "p.R399Q", "p.V262L", "p.V624I")
first_vec_vaf_063 <- c(49.56, 40.15, 37.7, 39.79, 44.91, 24.68, 67.77, 42.67, 68.13, 16.61)
first_vec_tmr_063 <- rep(0.88, length(first_vec_gene_063))
table_063_f <- as.data.frame(cbind(first_vec_sample_063, first_vec_gene_063, first_vec_base_change_063, first_vec_aa_change_063, first_vec_vaf_063, first_vec_tmr_063))

# Making one big list of tables for all of the genetic mutations, one table per patient
thetables <- list(table_2M_f, table_2N, table_3K, table_3K_f, table_4L_f, table_238, table_499, table_499_f, table_602, table_063_f)

# Assigning column names for each table in the list
for(thetable in thetables){
  colnames(thetable) <- c("Sample", "Gene", "Base_change", "AA_change", "VAF", "TMR")
  tmp <- thetable
  assign(paste0("table_", thetable[1,1]), tmp)
}

# Combining all of the tables into one big table
processed_NGS_table <- rbind(table_2M_f, table_2N, table_3K, table_3K_f, table_4L_f, table_pt238, table_pt499, table_pt499_f, table_pt602, table_pt063_f)

# Importing the non-batch corrected GeoMx data and metadata
vsd <- read.csv("outputs/preprocessing/vsd_not_batch_corrected_vsd.csv", row.names = 1)
metadata <- read.csv("outputs/preprocessing/metadata.csv", row.names = 1)

# 3K and 499 were sequenced twice. Removing duplicate mutations from the two different runs.
which_3k_duplicates <- processed_NGS_table %>% filter(Sample %in% c("3K_f", "3K")) %>% dplyr::select("AA_change") %>% as.vector() %>% unlist() %>% as.vector() %>% duplicated() %>% which() 
which_499_duplicates <- processed_NGS_table %>% filter(Sample %in% c("pt499_f", "pt499")) %>% dplyr::select("AA_change") %>% as.vector() %>% unlist() %>% as.vector() %>% duplicated() %>% which() 
total_3k <- processed_NGS_table %>% filter(Sample %in% c("3K_f", "3K"))
total_499 <- processed_NGS_table %>% filter(Sample %in% c("pt499_f", "pt499"))
unique_3K <- total_3k[-which_3k_duplicates, ]
unique_499 <- total_499[-which_499_duplicates, ]

# Making a single table with all of the NGS results from all of the patients
processed_unique_NGS_table <- rbind(table_2M_f, table_2N, unique_3K, table_4L_f, table_pt238, unique_499, table_pt602, table_pt063_f)
processed_unique_NGS_table$Gene_AA <- paste0(processed_unique_NGS_table$Gene, processed_unique_NGS_table$AA_change)

# Assigning patient aliases to the lab sample IDs for each tissue
processed_unique_NGS_table$Patient.Alias <- "tmp"
processed_unique_NGS_table$Patient.Alias[which(processed_unique_NGS_table$Sample == "2M_f")] <- 1
processed_unique_NGS_table$Patient.Alias[which(processed_unique_NGS_table$Sample == "2N")] <- 2
processed_unique_NGS_table$Patient.Alias[which(processed_unique_NGS_table$Sample == "3K")] <- 3
processed_unique_NGS_table$Patient.Alias[which(processed_unique_NGS_table$Sample == "4L_f")] <- 4
processed_unique_NGS_table$Patient.Alias[which(processed_unique_NGS_table$Sample == "pt238")] <- 5
processed_unique_NGS_table$Patient.Alias[which(processed_unique_NGS_table$Sample == "pt499")] <- 6
processed_unique_NGS_table$Patient.Alias[which(processed_unique_NGS_table$Sample == "pt499_f")] <- 6
processed_unique_NGS_table$Patient.Alias[which(processed_unique_NGS_table$Sample == "pt602")] <- 7
processed_unique_NGS_table$Patient.Alias[which(processed_unique_NGS_table$Sample == "pt063_f")] <- 8

# Assigning p53 status to each patient
processed_unique_NGS_table$stain_p53 <- "lost"
processed_unique_NGS_table$stain_p53[which(processed_unique_NGS_table$Patient.Alias %in% c(1,6))] <- "positive"
processed_unique_NGS_table$stain_p53[which(processed_unique_NGS_table$Patient.Alias %in% c(2,3,4,7))] <- "overexpressed"

# Assigning SMAD4 status to each patient
processed_unique_NGS_table$stain_SMAD4 <- "lost"
processed_unique_NGS_table$stain_SMAD4[which(processed_unique_NGS_table$Patient.Alias %in% c(2,4,8,5))] <- "positive"

# saving the organized NGS table
saveRDS(processed_unique_NGS_table, "inputs/combined_unique_NGS_table.rds")
write.csv(processed_unique_NGS_table, "inputs/combined_unique_NGS_table.csv")

metadata <- readRDS("outputs/preprocessing/metadata_with_VI_annotations.rds")

rownames(processed_unique_NGS_table) <- seq(1:nrow(processed_unique_NGS_table))
descending_order <- names(sort(table(processed_unique_NGS_table$Gene), decreasing = T))
stacked_patients <- processed_unique_NGS_table[, which(colnames(processed_unique_NGS_table) %in% c("Gene", "Patient.Alias"))]
stacked_patients$Gene <- factor(stacked_patients$Gene, levels = descending_order)

# FIGURE S9A ---- stacked bar plot for the four cancer genes of interest by patient alias ####
stacked_patients_filtered <- stacked_patients %>% filter(Gene %in% c("KRAS", "SMAD4", "CDKN2A", "TP53"))

pdf("outputs/plots_tables_objects/mutation_bar_plot_only_4_genes.pdf", width = 4, height = 6)
ggplot(stacked_patients_filtered, aes(x = Gene, fill = Patient.Alias)) +
  geom_bar(position = "stack") +
  labs(title = "Gene Mutations per Patient",
       x = "Gene Mutation",
       y = "Number of Patients",
       fill = "Patient Alias") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black", face = "italic"))
dev.off()
#####

# making stacked bar plot where the Y axis is the type of VI, X has four annotations. first row is p53 stain, second row is P53 mutation, third row is SMAD4 stain, fourth row is SMAD4 mutation

countlist <- list()
for(x in 1:8){
  VI_IN_Like_count <- metadata %>% filter(Patient.Alias == x & VI_Classification == "VI-IN-Like") %>% nrow()
  VI_Destructive_count <-  metadata %>% filter(Patient.Alias == x & VI_Classification == "VI-Destructive") %>% nrow()
  VI_Conventional_count <-  metadata %>% filter(Patient.Alias == x & VI_Classification == "VI-Conventional") %>% nrow()
  countlist[[x]] <- c(VI_IN_Like_count, VI_Destructive_count, VI_Conventional_count)
}

# making the table of VI type counts per patient
counttable <- as.data.frame(do.call(rbind, countlist))
colnames(counttable) <- c("VI-IN-Like", "VI-Destructive", "VI-Conventional")
rownames(counttable) <- seq(1:8)
counttable$Patient.Alias <- factor(c(1,2,3,4,5,6,7,8))

# Adding the SMAD4 and p53 staining data to the GeoMx metadata
alias_smad4_p53 <- unique(processed_unique_NGS_table[, c("Patient.Alias", "stain_p53", "stain_SMAD4")])
metadata <- merge(metadata, alias_smad4_p53, by = "Patient.Alias")

# Adding the mutational status for KRAS, CKDN2A, TP53, and SMAD4 to the GeoMx metadata
tmp <- unique(processed_unique_NGS_table[, c("Gene", "Base_change", "Patient.Alias", "AA_change")]) %>% filter(Gene %in% c("KRAS", "CDKN2A", "TP53", "SMAD4"))
tmp$mutation_category <- c("Missense", "Frameshift", "Missense", "Missense", "Missense", "Missense", "Missense", "Missense", "Missense", "Missense", "Missense", "Missense", "Missense", "Missense", "Missense", "Missense", "Missense", "Missense", "Nonsense")
tmp <- pivot_wider(tmp, names_from = "Gene", values_from = "mutation_category", values_fill = NA)
tmp <- as.data.frame(tmp)
tmp <- tmp[, c("Patient.Alias", "KRAS", "SMAD4", "TP53", "CDKN2A")]
tmp <- tmp %>%
  group_by(Patient.Alias) %>%
  summarise_at(vars(KRAS:CDKN2A), ~paste(na.omit(.), collapse = "; ")) %>%
  ungroup()
tmp <- as.data.frame(tmp)
tmp[tmp == ""] <- "WT"
colnames(tmp) <- c("Patient.Alias", "KRAS_mutation", "SMAD4_mutation", "TP53_mutation", "CDKN2A_mutation")
metadata <- merge(metadata, tmp, by = "Patient.Alias")

# getting a vector of the smad4 and p53 staining for each patient
tmp <- unique(metadata[, c("Patient.Alias","stain_p53", "stain_SMAD4")])
p53stain <- tmp$stain_p53
smad4stain <- tmp$stain_SMAD4

# making vectors for smad4 and p53 mutations
tmp <- unique(metadata[, c("Patient.Alias","SMAD4_mutation", "TP53_mutation")])
p53muts <- tmp$TP53_mutation
SMAD4muts <- tmp$SMAD4_mutation

annotations <- data.frame(
  Patient.Alias = rownames(counttable),
  p53Staining = p53stain,
  SMAD4Staining = smad4stain,
  p53Mutation = p53muts,
  SMAD4Mutation = SMAD4muts
)

colnames(annotations) <- c("Patient.Alias", "p53_Stain", "SMAD4_Stain", "p53_Mutation", "SMAD4_Mutation")

# Create the bar plot
bar_plot <- ggplot(melt(counttable, id.vars = "Patient.Alias"), aes(x = Patient.Alias, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of VI Types per patient", x = "Patient Alias", y = "Number of VI foci", fill = "Sample Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), title = element_text(hjust = 0.5)) + 
  scale_fill_manual(name = "Sample Type", labels = c("VI-IN-Like", "VI-Destructive", "VI-Conventional"), values = c("#36DAFF", "#9C27B0", "#A9A9A9"))

# Create annotations for p53 and SMAD4 staining and mutation levels
ha <- HeatmapAnnotation(
  `p53 Stain` = annotations$p53_Stain,
  `p53 Mutation` = annotations$p53_Mutation,
  `SMAD4 Stain` = annotations$SMAD4_Stain,
  `SMAD4 Mutation` = annotations$SMAD4_Mutation,
  col = list(
    `p53 Stain` = c("lost" = "black", "positive" = "#77AC4B", "overexpressed" = "gold"),
    `p53 Mutation` = c("WT" = "#77AC4B", "Missense" = "brown", "Nonsense" = "red"),
    `SMAD4 Stain` = c("lost" = "black", "positive" = "#77AC4B"),
    `SMAD4 Mutation` = c("Frameshift" = "pink", "Missense" = "brown", "WT" = "#77AC4B")
  )
)

# # Create a blank heatmap with annotations
annotation_heatmap <- Heatmap(matrix(nrow = 0, ncol = 8), bottom_annotation = ha, show_heatmap_legend = T)

# FIGURE S9B ---- VI type bar plot annotated by SMAD4 and p53 status ####
pdf("outputs/plots_tables_objects/smad4_p53_heatmap_annotations_consolidated.pdf")
draw(annotation_heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

pdf("outputs/plots_tables_objects/VI_type_bar_plot.pdf")
bar_plot
dev.off()
#####

# saving complete metadata table
saveRDS(metadata, "outputs/preprocessing/metadata_with_VI_subtypes_and_NGS.rds")
write.csv(metadata, "outputs/preprocessing/metadata_with_VI_subtypes_and_NGS.csv")
rownames(metadata) <- metadata$Name
metadata$Type <- factor(metadata$Type, levels=c('ND','PDAC', 'PNI','VI'))


# Completing VI vs PDAC in SMAD4-WT and SMAD4 Mutant samples

# importing non batch corrected GeoMx data
noncorrected_vsd <- read.csv("outputs/preprocessing/vsd_not_batch_corrected_vsd.csv", row.names = 1)

# making pheno data and expression set for limma
vsd <- as.matrix(vsd)
pheno = new("AnnotatedDataFrame", data=metadata)
eset = ExpressionSet(assayData = vsd, phenoData = pheno)
edata = Biobase::exprs(eset)

# batch correction on batch
batch_adjusted <- sva::ComBat(dat = edata, batch = metadata$Batch, par.prior = TRUE)

# labeling smad4 as normal/binary if either mutated or lost/overexpressed
metadata$smad4_total_binary <- "normal"
metadata$smad4_total_binary[which(metadata$SMAD4_mutation != "WT" | metadata$stain_SMAD4 != "positive")] <- "aberrant"

# separating the metadata by SMAD4 groups
smad4_normal_metadata <- metadata %>% filter(smad4_total_binary == "normal")
smad4_aberrant_metadata <- metadata %>% filter(smad4_total_binary == "aberrant")
smad4_normal_metadata$Type <- factor(smad4_normal_metadata$Type, levels=c('ND','PDAC', 'PNI','VI'))
smad4_aberrant_metadata$Type <- factor(smad4_aberrant_metadata$Type, levels=c('ND','PDAC', 'PNI','VI'))

# separating the GeoMx data by SMAD4 groups
batch_adjusted_smad4_aberrant <- batch_adjusted[,which(metadata$smad4_total_binary == "aberrant")]
batch_adjusted_smad4_normal <- batch_adjusted[,which(metadata$smad4_total_binary == "normal")]

# batch correction on patient for each data set broken up by SMAD4 group
batch_adjusted_smad4_aberrant <- sva::ComBat(dat = batch_adjusted_smad4_aberrant, 
                                             batch =smad4_aberrant_metadata$Slide, par.prior = TRUE)

batch_adjusted_smad4_normal <- sva::ComBat(dat = batch_adjusted_smad4_normal, 
                                           batch =smad4_normal_metadata$Slide, par.prior = TRUE)

# making pheno data and expression set for limma for SMAD4 normal
vsd_normal <- as.matrix(batch_adjusted_smad4_normal)
pheno = new("AnnotatedDataFrame", data=smad4_normal_metadata)
eset = Biobase::ExpressionSet(assayData = vsd_normal, phenoData = pheno)
edata = Biobase::exprs(eset)

# Making matrix design
design = model.matrix(~0 + pheno$Type)
colnames(design) <- levels(pheno$Type)

# Making contrasts for differential expression
contr.matrix <- makeContrasts(
  VIvsPDAC = VI - PDAC,
  levels = colnames(design))

# doing limma differential expression on the batch adjusted matrix
eset = Biobase::ExpressionSet(assayData = vsd_normal, phenoData = pheno)
fit <- lmFit(eset, design)
vfit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(vfit)

vi.vs.pdac <- topTreat(efit, coef=1, n=Inf)
volcanoPlot(vi.vs.pdac, title_1 = "VI", title_2 = "PDAC", batch = "smad4_normal", directory = "outputs/plots_tables_objects/smad4normal/")

# making pheno data and expression set for limma
vsd_aberrant <- as.matrix(batch_adjusted_smad4_aberrant)
pheno = new("AnnotatedDataFrame", data=smad4_aberrant_metadata)
eset = Biobase::ExpressionSet(assayData = vsd_aberrant, phenoData = pheno)
edata = Biobase::exprs(eset)
pheno$Type
# Making matrix design
design = model.matrix(~0 + pheno$Type)
colnames(design) <- levels(pheno$Type)

# Making contrasts for differential expression
contr.matrix <- makeContrasts(
  VIvsPDAC = VI - PDAC,
  levels = colnames(design))

# doing limma differential expression on the batch adjusted matrix
eset = Biobase::ExpressionSet(assayData = vsd_aberrant, phenoData = pheno)
fit <- lmFit(eset, design)
vfit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(vfit)

vi.vs.pdac <- topTreat(efit, coef=1, n=Inf)
volcanoPlot(vi.vs.pdac, title_1 = "VI", title_2 = "PDAC", batch = "smad4_aberrant", directory = "outputs/plots_tables_objects/smad4aberrant/")

# Making the venn diagrams
smad4_normal_vi_vs_pdac <- read.csv("outputs/plots_tables_objects/smad4normal/smad4_normal_VIvsPDAC_all_genes.csv")
smad4_aberrant_vi_vs_pdac <- read.csv("outputs/plots_tables_objects/smad4aberrant/smad4_aberrant_VIvsPDAC_all_genes.csv")

# setting row names
rownames(smad4_normal_vi_vs_pdac) <- smad4_normal_vi_vs_pdac$rn
rownames(smad4_aberrant_vi_vs_pdac) <- smad4_aberrant_vi_vs_pdac$rn

# setting equivalent row order between smad4 aberrant and normals
smad4_normal_vi_vs_pdac <- smad4_normal_vi_vs_pdac[sort(smad4_normal_vi_vs_pdac$rn, decreasing = F),]
smad4_aberrant_vi_vs_pdac <- smad4_aberrant_vi_vs_pdac[sort(smad4_normal_vi_vs_pdac$rn, decreasing = F),]

# subtracting the log fold change of SMAD4 normal from SMAD4 aberrant
the_lfc_difference <- smad4_aberrant_vi_vs_pdac$logFC - smad4_normal_vi_vs_pdac$logFC

# naming the vector
names(the_lfc_difference) <- smad4_normal_vi_vs_pdac$rn

# sorting the genes by the magnitude of the difference between LFC in the two groups
the_diff <- sort(abs(the_lfc_difference), decreasing = T)

# Getting all genes which have greater than one log fold change difference between SMAD4-aberrant and SMAD4-normal
the_diffs <- which(the_diff > 1)

# Getting genes increased in VI in both groups
normal_vi_genes <- smad4_normal_vi_vs_pdac %>% filter(Group == "Uppers")
aberrant_vi_genes <- smad4_aberrant_vi_vs_pdac %>% filter(Group == "Uppers")

# Getting genes increased in PDAC in both groups
normal_pdac_genes <- smad4_normal_vi_vs_pdac %>% filter(Group == "Lowers")
aberrant_pdac_genes <- smad4_aberrant_vi_vs_pdac %>% filter(Group == "Lowers")

# SMAD4_WT VI genes
SMAD4_WT_VI_Markers <- setdiff(normal_vi_genes$rn, aberrant_vi_genes$rn)[setdiff(normal_vi_genes$rn, aberrant_vi_genes$rn) %in% names(the_diffs)]

# SMAD4 Aberrant VI genes
SMAD4_Aberrant_VI_Markers <- setdiff(aberrant_vi_genes$rn, normal_vi_genes$rn)[setdiff(aberrant_vi_genes$rn, normal_vi_genes$rn) %in% names(the_diffs)]

# SMAD4 WT PDAC Genes
SMAD4_WT_PDAC_Markers <- setdiff(normal_pdac_genes$rn, aberrant_pdac_genes$rn)[setdiff(normal_pdac_genes$rn, aberrant_pdac_genes$rn) %in% names(the_diffs)]

# SMAD4 Aberrant PDAC genes
SMAD4_Aberrant_PDAC_Markers <- setdiff(aberrant_pdac_genes$rn, normal_pdac_genes$rn)[setdiff(aberrant_pdac_genes$rn, normal_pdac_genes$rn) %in% names(the_diffs)]

# FIGURE S9D ---- venn diagrams for the VI-PDAC DEGs in both SMAD4 groups ####
# Making venn diagrams for the overlaps
venn.diagram(x = list(normal_vi_genes$rn, aberrant_vi_genes$rn), filename = "outputs/plots_tables_objects/smad4_vi_aberrant_normal.png", category.names = c("", ""),
             fill = c("#D6A857", "#2196F3"), main = "VI enriched")

venn.diagram(x = list(normal_pdac_genes$rn, aberrant_pdac_genes$rn), filename = "outputs/plots_tables_objects/pdac_aberrant_normal.png", category.names = c("", ""),
             fill = c("#D6A857", "#2196F3"), main = "PDAC enriched")

venn.diagram(x = list(SMAD4_WT_VI_Markers, SMAD4_Aberrant_VI_Markers), filename = "outputs/plots_tables_objects/vi_aberrant_normal_post_filter.png", category.names = c("", ""),
             fill = c("#D6A857", "#2196F3"), main = "VI_enriched")

venn.diagram(x = list(SMAD4_WT_PDAC_Markers, SMAD4_Aberrant_PDAC_Markers), filename = "outputs/plots_tables_objects/pdac_aberrant_normal_post_filter.png", category.names = c("", ""),
             fill = c("#D6A857", "#2196F3"), main = "PDAC_enriched")
#####

# All of the categories for the gene groups
namelist <- list("wt_vi", "ab_vi", "wt_pdac", "ab_pdac")
genelist <- list(SMAD4_WT_VI_Markers, SMAD4_Aberrant_VI_Markers, SMAD4_WT_PDAC_Markers, SMAD4_Aberrant_PDAC_Markers)

# making table of p values for each comparison
vipdac_normal_p_table <- data.frame("Genes" = smad4_normal_vi_vs_pdac$rn, "group_1" = "VI", "group_2" = "PDAC", "p" = smad4_normal_vi_vs_pdac$adj.P.Val, "LFC" = smad4_normal_vi_vs_pdac$logFC, "SMAD4" = "normal")
vipdac_aberrant_p_table <- data.frame("Genes" = smad4_aberrant_vi_vs_pdac$rn, "group_1" = "VI", "group_2" = "PDAC", "p" = smad4_aberrant_vi_vs_pdac$adj.P.Val, "LFC" = smad4_aberrant_vi_vs_pdac$logFC, "SMAD4" = "aberrant")

# Reading in the fully batch corrected (by batch and by patient)
vsd <- read.csv("outputs/preprocessing/fully_batch_corrected_vsd.csv", row.names = 1)
metadata$smad4_total_binary <- factor(metadata$smad4_total_binary, levels = c("normal", "aberrant"))

# FIGURE S9E --- getting both plots of the VI and PDAC marker genes for both SMAD4 groups####
p_table <- rbind(vipdac_normal_p_table, vipdac_aberrant_p_table)
for(x in 1:4){
  genesSig <- genelist[[x]]
  metadataDat <- data.frame(metadata,t(vsd[genesSig,]))
  metadataDat <- metadataDat %>% filter(Type == "PDAC" | Type == "VI")
  for (g in genesSig) {
    metadataDatgene <- metadataDat %>% dplyr::select("Type", g)
    data_for_p_val_manual <- p_table %>% filter(Genes == make.names(g)) %>% dplyr::select("group_1", "group_2", "p", "LFC", "SMAD4")
    data_for_p_val_manual <- data_for_p_val_manual %>% dplyr::select("group_1", "group_2", "p", "SMAD4")
    colnames(data_for_p_val_manual) <- c("group1", "group2", "p.adj", "smad4_total_binary")
    datman <- tibble(data_for_p_val_manual)
    
    # Convert p-values to asterisks
    datman$p.signif <- sapply(datman$p.adj, p_to_asterisk)
    
    yposition <- max(metadataDatgene[,make.names(g)]) * 1.05
    stat.test <- stat_pvalue_manual(datman, x = 'smad4_total_binary', y.position = yposition, label = "p.signif")
    pdf(paste0("outputs/geneplots/smad4/",namelist[x],"/", g,"_boxplot.pdf"), width = 4, height = 6)
    print(ggplot(metadataDat, aes_string(x='smad4_total_binary', y=make.names(g), color='Type')) +
            geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
            labs(color = "SMAD4 status", x = "", y = "Normalized expression", title = g) +
            scale_color_manual(name = "Tissue type", labels = c("PDAC", "VI"), values = c("#DF2E2E","#2196F3")) +
            theme(axis.text.x = element_text(color = "black",
                                             size = 12, angle = 45, hjust = 1),
                  plot.title = element_text(hjust = 0.5, size = 18, face = "italic")) +
            stat.test)
    dev.off()
  }
}
#####

# Doing GSEA

# doing regular fgsea with the custom genesets

# reading in the cell types from Hwang et al
cell_types <- read.csv("inputs/snRNA_seq_cell_type_sets.csv")
cell_states <- read.csv("inputs/snRNA_seq_cell_state_sets.csv")
cell_types_list <- lapply(cell_types, as.vector)
cell_states_list <- lapply(cell_states, as.vector)

# Reading in the cogaps pattern markers
jp4thresh <- readRDS('outputs/thresholded_patterns.rds')

# Separating the pattern markers by pattern
Pattern_1 <- jp4thresh$PatternMarkers[[1]]
Pattern_2 <- jp4thresh$PatternMarkers[[2]]
Pattern_3 <- jp4thresh$PatternMarkers[[3]]
Pattern_4 <- jp4thresh$PatternMarkers[[4]]

# The invasion types and Moffitt types
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

# Putting all of the gene sets into a single named list
genesets <- list(basal, classical, collective, mesenchymal, Pattern_1, Pattern_3, Pattern_4, cell_types_list$Acinar.like, cell_types_list$Basaloid, cell_types_list$Squamoid, cell_types_list$Mesenchymal, cell_types_list$Neuroendocrine.like, cell_types_list$Neural.like.progenitor, cell_states_list$Cycling..S., cell_states_list$Cycling..G2.M.,
                 cell_states_list$MYC.signaling, cell_states_list$Adhesive, cell_states_list$Ribosomal, cell_states_list$Interferon.signaling, cell_states_list$TNF.NFkB.signaling)
names(genesets) <- c("moffitt_basal", "moffitt_classical", "wood_collective", "wood_mesenchymal","Pattern_1", "Pattern_3", "Pattern_4", "acinar.like", "basaloid", "squamoid", "mesenchymal", "neuroendocrine.like", "neural.like.progenitor", 
                     "cycling.s", "cycling.G2M", "MYC.signaling", "adhesive", "ribosomal", "interferon", "TNF.NFkB.signaling")

# the LFCs of VI vs PDAC DEGs in SMAD4-WT
normal_fgsea <- smad4_normal_vi_vs_pdac$logFC
names(normal_fgsea) <- smad4_normal_vi_vs_pdac$rn

# running fgsea on these genes using LFC as the ordering statistic
tmp_normalfgsea <- fgsea(pathways = genesets, stats = normal_fgsea)

# turning results into data frame
tmp_normalfgsea <- as.data.frame(tmp_normalfgsea)

# Ordering genes by NES
NES_normal <- tmp_normalfgsea[order(tmp_normalfgsea$NES, decreasing = T),]

# Performing the same steps in SMAD4-aberrant samples
aberrant_fgsea <- smad4_aberrant_vi_vs_pdac$logFC
names(aberrant_fgsea) <- smad4_aberrant_vi_vs_pdac$rn
tmp_aberrantfgsea <- fgsea(pathways = genesets, stats = aberrant_fgsea)
tmp_aberrant_fgsea <- as.data.frame(tmp_aberrantfgsea)
NES_aberrant <- tmp_aberrantfgsea[order(tmp_aberrantfgsea$NES, decreasing = T),]

# making lollipop plot
positives <- NES_normal %>% filter(NES > 0)
negatives <- NES_normal %>% filter(NES < 0)

# making lollipop plot
positives <- positives[order(positives$NES, decreasing = F), ]
negatives <- negatives[order(negatives$NES, decreasing = F), ]

# makign lollipop plot
positives <- positives %>% mutate(group = "positives")
negatives <- negatives %>% mutate(group = "negatives")

# making lollipop plot
thedata <- rbind(negatives, positives)
thedata <- thedata %>% dplyr::select("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "group")
thedata$pathway <- factor(thedata$pathway, levels = thedata$pathway)

# FIGURE S9C ---- pathway analysis lollipop plot for SMAD4-normal VI-PDAC comparison ####
pdf("outputs/plots_tables_objects/patterns_SMAD4_normal_lollipop.pdf", width = 4, height = 4)
ggplot(thedata, aes(x=pathway, y=NES)) +
  geom_segment(aes(x=pathway, xend=pathway, y=0, yend=NES), size=0.5, alpha=0.5) +
  geom_point(size = 3, aes(fill= group), alpha=1, shape=21, stroke=1.4) +
  scale_fill_manual(values = c("#DF2E2E", "#2196F3")) +
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

# making lollipop plot
positives <- NES_aberrant %>% filter(NES > 0)
negatives <- NES_aberrant %>% filter(NES < 0)

# making lollipop plot
positives <- positives[order(positives$NES, decreasing = F), ]
negatives <- negatives[order(negatives$NES, decreasing = F), ]

# makign lollipop plot
positives <- positives %>% mutate(group = "positives")
negatives <- negatives %>% mutate(group = "negatives")

# making lollipop plot
thedata <- rbind(negatives, positives)
thedata <- thedata %>% dplyr::select("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "group")
thedata$pathway <- factor(thedata$pathway, levels = thedata$pathway)

# FIGURE S9C ---- pathway analysis lollipop plot for SMAD4-aberrant VI-PDAC comparison ####
pdf("outputs/plots_tables_objects/patterns_SMAD4_aberrant_lollipop.pdf", height = 4, width = 4)
ggplot(thedata, aes(x=pathway, y=NES)) +
  geom_segment(aes(x=pathway, xend=pathway, y=0, yend=NES), size=0.5, alpha=0.5) +
  geom_point(size = 3, aes(fill= group), alpha=1, shape=21, stroke=1.4) +
  scale_fill_manual(values = c("#DF2E2E", "#2196F3")) +
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

