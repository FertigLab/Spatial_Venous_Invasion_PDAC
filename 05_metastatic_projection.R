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
library(epitools)
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

# this data was created in the "left_over_icgcs" script
left_overs <- read.csv("outputs/left_over_icgc_samples.csv")
rownames(left_overs) <- left_overs$X
left_overs <- left_overs[,-1]
mergeddat <- merge(as.data.frame(thedata), as.data.frame(left_overs), by = 'row.names')
rownames(mergeddat) <- mergeddat[,1]
mergeddat <- mergeddat[,-1]

# filtering for the data present in paired_unpaired (from the paper)
paired_unpaired <- read.csv("inputs/paired_unpaired.csv", header = T)

# making a vector of all the IDs from the paper 
paper_idents <- paired_unpaired$RNASeq.id[!is.na(paired_unpaired$RNASeq.id)]

# identifying typos in the paper idents table
typos <- paper_idents[which(!paper_idents %in% colnames(mergeddat))]
fixed <- c("PCSI_0632_Lv_M_526", "PCSI_0663_Lv_M_526", "PCSI_0664_Lv_M_526")

# fixing the typos in the paired-unpaired data table
paired_unpaired$RNASeq.id[which(paired_unpaired$RNASeq.id %in% typos)] <- fixed

# eliminating the samples without RNA Seq from the data set & xenograft and low cellularity
paired_unpaired <- paired_unpaired %>% filter(!is.na(RNASeq.id) & tissue.source != "xenograft, fresh frozen" & cellularity.celluloid >= 0.70)

# getting updated/filtered vector of sample IDs
paper_idents <- paired_unpaired$RNASeq.id

# filtering out all of the sample IDs not present in the paper's list
thedata <- mergeddat %>% dplyr::select(all_of(paper_idents))

#write.csv(paired_unpaired, inputs/paired_unpaired_notypos.csv")

# turning ENSMBL into HGNC
thedata <- read.csv("inputs/raw_matrix_all_from_paper.csv", row.names = 1)
paired_unpaired <- read.csv("inputs/paired_unpaired_notypos.csv")

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

# making my own phenodata
pheno <- data.frame()
patient_id<- str_before_nth(colnames(newdata), "_", 2)
primemets <- rep("primary",times=ncol(newdata))
primemets[which(!is.na(str_match(colnames(newdata), "_M_")))] <- "metastatic"
primemets[which(!is.na(str_match(colnames(newdata), "_M")))] <- "metastatic"
organ <- rep("liver",times=ncol(newdata))
organ[which(!is.na(str_match(colnames(newdata), "Lu")))] <- "lung"
organ[which(!is.na(str_match(colnames(newdata), "Ab")))] <- "abdominal"
organ[which(!is.na(str_match(colnames(newdata), "Ag")))] <- "adrenal_gland"
organ[which(!is.na(str_match(colnames(newdata), "Di")))] <- "diaphragm"
organ[which(!is.na(str_match(colnames(newdata), "Du")))] <- "duodenum"
organ[which(!is.na(str_match(colnames(newdata), "Hr")))] <- "heart"
organ[which(!is.na(str_match(colnames(newdata), "Ln")))] <- "lymph_node"
organ[which(!is.na(str_match(colnames(newdata), "Mu")))] <- "muscle"
organ[which(!is.na(str_match(colnames(newdata), "Om")))] <- "omentum"
organ[which(!is.na(str_match(colnames(newdata), "Pa")))] <- "pancreas"
organ[which(!is.na(str_match(colnames(newdata), "Si")))] <- "small_intestine"
organ[which(!is.na(str_match(colnames(newdata), "Sp")))] <- "spleen"
organ[which(!is.na(str_match(colnames(newdata), "St")))] <- "stomach"
organ[which(!is.na(str_match(colnames(newdata), "Pm")))] <- "peritoneum"

# making metadata
mypheno <- as.data.frame(cbind(patient_id, primemets, organ))
colnames(mypheno) <- c("patient_id", "stage", "organ")
rownames(mypheno) <- colnames(newdata)
rownames(paired_unpaired) <- paired_unpaired$RNASeq.id
mergedpheno <- merge(paired_unpaired, mypheno, by = 0)
mergedpheno <- mergedpheno %>% filter(organ != "lymph_node")

# cleaning the data matrix to include same samples as metadata
newdata <- newdata[,which(colnames(newdata) %in% mergedpheno$RNASeq.id)]

prims <- mergedpheno %>% filter(stage == "primary")
mets <- mergedpheno %>% filter(stage == "metastatic")

# removing low count genes from matrix
table(rowMeans(newdata))
myCPM <- edgeR::cpm(newdata)

thresh <- myCPM > 1
table(rowSums(thresh))

keep <- rowSums(thresh) >= ncol(newdata)/10
dim(newdata)
newdata <- newdata[keep,]
cogaps <- readRDS("inputs/batch_merged_cogapsP4.rds")

# keeping row names from geomx
geomxgenes <- rownames(cogaps@loadingStdDev)
newdata <- na.omit(newdata[geomxgenes,])


## Normalizing with VST
dds <- DESeqDataSetFromMatrix(countData = newdata,
                              colData = mergedpheno,
                              design= ~0 + stage)
vsd <- assay(vst(dds, blind=F))
vsd <- as.data.frame(vsd)

pmarkers <- patternMarkers(cogaps, threshold="cut")
project_embed <- projectR(data = as.matrix(vsd), loadings = cogaps)

# making table of the pattern weights
pattern_weights <- as.data.frame(t(project_embed))

# gathering a table of patterns and the score
xyz <- gather(pattern_weights, key = "Pattern", value = "score")

# making a vector of the patient IDs, repeating 4x in the order that they were originally listed in the metadata
the_idents <- c(rownames(pattern_weights), rownames(pattern_weights), rownames(pattern_weights), rownames(pattern_weights))

# making a column in pattern weights of the patient IDs
pattern_weights$identities <- rownames(pattern_weights)

# adding the repeating patient IDs to the gathered table
the_table <- cbind(xyz, the_idents)

# adding the patient IDs to the metadata as a column
mergedpheno$identities <- mergedpheno$Row.names

# merging the metadata with the pattern weights by the patient identities
themetadata <- merge(mergedpheno, pattern_weights, by = "identities")

# making the metadata rownames the patient identities
rownames(themetadata) <- themetadata$identities

# renaming the patient identities
colnames(the_table)[3] <- "identities"

# joining the gathered table with metadata to boxplot it
final_table <- left_join(the_table, themetadata, by='identities')

# only filtering to Pattern 1 because I don't want duplicate samples (doesn't necessarily need to be pattern 1). I'm going to use this table to compare the patterns between mets and primary.
prims_for_wilcox <- final_table %>% filter(stage == "primary" & Pattern == "Pattern_1")
mets_for_wilcox <- final_table %>% filter(stage == "metastatic" & Pattern == "Pattern_1")

# doing stats for the CoGAPS pattern projection
metprim1 <- wilcox.test(prims_for_wilcox$Pattern_1, mets_for_wilcox$Pattern_1)
metprim1 <- metprim1$p.value
metprim1 <- data.frame("Pattern" = "1", "group1" = "metastatic", "group2" = "primary", "p" = metprim1)
metprim2 <- wilcox.test(prims_for_wilcox$Pattern_2, mets_for_wilcox$Pattern_2)
metprim2 <- metprim2$p.value
metprim2 <- data.frame("Pattern" = "2", "group1" = "metastatic", "group2" = "primary", "p" = metprim2)
metprim3 <- wilcox.test(prims_for_wilcox$Pattern_3, mets_for_wilcox$Pattern_3)
metprim3 <- metprim3$p.value
metprim3 <- data.frame("Pattern" = "3", "group1" = "metastatic", "group2" = "primary", "p" = metprim3)
metprim4 <- wilcox.test(prims_for_wilcox$Pattern_4, mets_for_wilcox$Pattern_4)
metprim4 <- metprim4$p.value
metprim4 <- data.frame("Pattern" = "4", "group1" = "metastatic", "group2" = "primary", "p" = metprim4)

# combining all the p values together in a single table
p_table <- rbind(metprim1, metprim2, metprim3, metprim4)

# FIGURE 7A ---- CoGAPS Pattern projection boxplots ####
for(x in 1:length(rownames(project_embed))){
  pdf(paste0("outputs/plots_tables_objects/pattern_",x,".pdf"), width = 4, height = 6)
  filtered_table <- final_table %>% filter(Pattern == rownames(project_embed)[x])
  data_for_p_val_manual <- p_table %>% filter(Pattern == x) %>% dplyr::select("group1", "group2", "p")
  datman <- tibble(data_for_p_val_manual)
  
  # Convert p-values to asterisks
  datman$p.signif <- sapply(datman$p, p_to_asterisk)
  
  datman$p <- signif(datman$p,digits=3)
  yposition <- max(filtered_table[,rownames(project_embed)[x]])*1.05
  stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
  print(ggplot(filtered_table, aes(x=stage, y=score, color = stage)) + 
          geom_boxplot() +
          geom_jitter(shape=16, position=position_jitter(0.2)) +
          scale_color_manual(name = "Stage", labels = c("Metastasis", "Primary"), values = c("#E81E63", "#3F51B5")) +
          scale_x_discrete(name = paste0("Pattern ",x), labels = c("", "")) +
          labs(title = paste0("Pattern ", x)) +
          theme(plot.title = element_text(hjust = 0.5, size = 18)) +
          stat.test) 
  dev.off()
}
#####

# Doing differential expression between metastasis and primary
rownames(mergedpheno) <- mergedpheno$Row.names
vsd <- vsd[,match(rownames(mergedpheno), colnames(vsd))]

# # Making matrix design
mergedpheno$stage <- factor(mergedpheno$stage, levels = c("primary", "metastatic"))
design = model.matrix(~0 + mergedpheno$stage)
colnames(design) <- c("primary", "metastatic")

# Making contrasts for differential expression comparisons
contr.matrix <- makeContrasts(
  metastatic_vs_primary = metastatic - primary,
  levels = colnames(design))

pheno = new("AnnotatedDataFrame", data=mergedpheno)
eset = ExpressionSet(assayData = as.matrix(vsd), phenoData = pheno)
fit <- lmFit(eset, design)
vfit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(vfit)
results_efit<-decideTests(efit, adj.P.val=0.05)
summary_efit <- summary(results_efit)

# comparing metastasis vs primary
met.vs.prim <- topTreat(efit, coef=1,, n=Inf)

# FIGURE S10A ---- volcano plot between metastasis and primary ####
volcanoPlot(met.vs.prim, title_1 = "Metastasis", title_2 = "Primary", batch = "vsd", upcolor ="#E81E63", downcolor = "#3F51B5", directory = "outputs/plots_tables_objects/")
#####

# reading in the metastasis vs primary DEGs
metprim <- read.csv("outputs/plots_tables_objects/vsd_MetastasisvsPrimary_DEGS.csv")

# reading in pattern 3 markers
VIgenes <- pmarkers$PatternMarkers[[3]]

# reading in pattern 1 markers
pdacgenes <- pmarkers$PatternMarkers[[1]]

# getting the pattern 3 markers which are also found in the metastasis vs primary DEGs
VIgenes <- VIgenes[which(VIgenes %in% rownames(met.vs.prim))]

# getting the pattern 1 markers which are also metastasis vs primary DEGs
pdacgenes <- pdacgenes[which(pdacgenes %in% rownames(met.vs.prim))]

# getting the genes elevated in metastasis
metgenes <- metprim %>% filter(Group == "Uppers")
metgenes <- metgenes$Genes

# getting the genes elevated in primary
primarygenes <- metprim %>% filter(Group == "Lowers")
primarygenes <- primarygenes$Genes

# intersecting the pattern 3 markers and the genes elevated in primary
VIprim <- intersect(VIgenes, primarygenes)

# intersecting the pattern 1 markers and the genes elevated in metastasis
pdacmet <- intersect(pdacgenes, metgenes)

# intersecting the pattern 3 markers and the genes elevated in metastasis
VImet <- intersect(VIgenes, metgenes)

# intersecting the pattern 1 markers and the genes elevated in priamary
pdacprim <- intersect(pdacgenes, primarygenes)

# Performing fisher's exact test for the bar plot

# getting number of genes in each category
pat1nonprim <- length(pdacgenes) - length(pdacprim)
pat1nonmet <- length(pdacgenes) - length(pdacmet)
pat1nonsig <- length(pdacgenes) - length(pdacprim) - length(pdacmet)

pat3nonprim <- length(VIgenes) - length(VIprim)
pat3nonmet <- length(VIgenes) - length(VImet)
pat3nonsig <- length(VIgenes) - length(VIprim) - length(VImet)

# making the outcomes table for the odds ratios/fisher's test
prim_outcomes <- c(rep(1, length(pdacprim)), rep(0, pat1nonprim), rep(1, length(VIprim)), rep(0, pat3nonprim))
met_outcomes <- c(rep(1, length(pdacmet)), rep(0, pat1nonmet), rep(1, length(VImet)), rep(0, pat3nonmet))

pat1_vs_pat3_primaries <- table(
  Population = c(rep("pat1", length(pdacgenes)), rep("pat3", length(VIgenes))),
  Outcome = prim_outcomes)

pat1_vs_pat3_mets <- table(
  Population = c(rep("pat1", length(pdacgenes)), rep("pat3", length(VIgenes))),
  Outcome = met_outcomes)

# the odds ratios and p values for the primary genes
theodds_prims <- oddsratio(pat1_vs_pat3_primaries, method = "fisher")

# the odds ratios and p values for the metastasis genes
theodds_mets <- oddsratio(pat1_vs_pat3_mets, method = "fisher")

the_df <- data.frame(matrix(NA, nrow = 256, ncol = 1))
the_df$labeled_geomx <- c(rep("Pattern_1", length(pdacgenes)), rep("Pattern_3", length(VIgenes)))
the_df$labeled_PACA <- c(rep("primary", length(pdacprim)), rep("neither", pat1nonsig), rep("metastasis", length(pdacmet)), rep("primary", length(VIprim)), rep("neither", pat3nonsig), rep("metastasis", length(VImet)))
the_df$count <- 1

# FIGURE 7B ---- making the bar plot for the pattern markers ####
pdf("outputs/plots_tables_objects/barplot_patternmarkers.pdf")
ggplot(the_df) + 
  geom_bar(aes(fill=labeled_PACA, y=count, x=labeled_geomx), position="stack", stat="identity") + 
  scale_fill_manual(name = "Primary-Metastasis RNA-Seq", labels = c("Significantly upregulated in metastasis", "Not signficantly upregulated", "Significantly upregulated in primary PDAC"), values = c("#E81E63", "#9C9C9C", "#3F51B5")) +
  scale_x_discrete(name = "GeoMx", labels = c("Pattern 1 Markers", "Pattern 3 Markers")) + 
  ylab("Number of Significantly Upregulated Genes")
dev.off()
#####

# making data frame to use to make box plots
data_for_plot <- as.data.frame(cbind(as.data.frame(t(vsd)), mergedpheno$stage, mergedpheno$organ))
colnames(data_for_plot)[which(colnames(data_for_plot) %in% c("mergedpheno$stage", "mergedpheno$organ"))] <- c("stage", "organ")
data_for_plot$stage <- factor(data_for_plot$stage, levels = c("metastatic", "primary"))
met_prim <- data.frame("Genes" = rownames(met.vs.prim), "group_1" = "metastatic", "group_2" = "primary", "p" = met.vs.prim$adj.P.Val)

# FIGURE 7C, S10B, S10D ---- boxplots of pattern markers which are also DEGs between mets and primaries ####
# making box plots of the pattern 1 markers which are elevated in metastasis relative to primary
for(gene in pdacmet){
  data_for_p_val_manual <- met_prim %>% filter(Genes == gene) %>% dplyr::select("group_1", "group_2", "p")
  colnames(data_for_p_val_manual) <- c("group1", "group2", "p.adj")
  datman <- tibble(data_for_p_val_manual)
  
  # Convert p-values to asterisks
  datman$p.signif <- sapply(datman$p.adj, p_to_asterisk)
  
  datman$p.adj <- signif(datman$p.adj,digits=3)
  yposition <- max(data_for_plot[,gene]) * 1.05
  stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
  pdf(paste0("outputs/geneplots/pdacmet/",gene,".pdf"), width = 4, height = 6)
  print(ggplot(data_for_plot, aes_string(x="stage", y=gene, color = "stage")) + 
          geom_boxplot() + 
          geom_jitter(shape=16, position=position_jitter(0.2)) +
          scale_x_discrete(labels = c("Metastasis", "Primary"), name = "Stage") +
          scale_color_manual(labels = c("Metastasis", "Primary"), name = "Stage", values = c("#E81E63", "#3F51B5"))+
          labs(y = "Normalized expression", title = gene) + 
          theme(plot.title = element_text(hjust = 0.5, size = 18, face = "italic")) + 
          stat.test)
  dev.off()
}

# making box plots of the pattern 3 markers which are elevated in primary relative to metastasis
for(gene in VIprim){
  data_for_p_val_manual <- met_prim %>% filter(Genes == gene) %>% dplyr::select("group_1", "group_2", "p")
  colnames(data_for_p_val_manual) <- c("group1", "group2", "p.adj")
  datman <- tibble(data_for_p_val_manual)
  
  # Convert p-values to asterisks
  datman$p.signif <- sapply(datman$p.adj, p_to_asterisk)
  
  datman$p.adj <- signif(datman$p.adj,digits=3)
  yposition <- max(data_for_plot[,gene]) * 1.05
  stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
  pdf(paste0("outputs/geneplots/viprim/",gene,".pdf"), width = 4, height = 6)
  print(ggplot(data_for_plot, aes_string(x="stage", y=gene, color = "stage")) + 
          geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))+  
          scale_x_discrete(labels = c("Metastasis", "Primary"), name = "Stage") +
          scale_color_manual(labels = c("Metastasis", "Primary"), name = "Stage", , values = c("#E81E63", "#3F51B5")) +
          labs(y = "Normalized expression", title = gene) + 
          theme(plot.title = element_text(hjust = 0.5, size = 18, face = "italic")) + 
          stat.test)
  dev.off()
}
#####
# FIGURE S10C ####
# making box plots of the pattern 1 markers which are elevated in primary relative to metastasis
for(gene in pdacprim){
  data_for_p_val_manual <- met_prim %>% filter(Genes == gene) %>% dplyr::select("group_1", "group_2", "p")
  colnames(data_for_p_val_manual) <- c("group1", "group2", "p.adj")
  datman <- tibble(data_for_p_val_manual)
  
  # Convert p-values to asterisks
  datman$p.signif <- sapply(datman$p.adj, p_to_asterisk)
  
  datman$p.adj <- signif(datman$p.adj,digits=3)
  yposition <- max(data_for_plot[,gene]) * 1.05
  stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
  pdf(paste0("outputs/geneplots/pdacprim/",gene,".pdf"), width = 4, height = 6)
  print(ggplot(data_for_plot, aes_string(x="stage", y=gene, color = "stage")) + 
          geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
          scale_x_discrete(labels = c("Metastasis", "Primary"), name = "Stage") +
          scale_color_manual(labels = c("Metastasis", "Primary"), name = "Stage", , values = c("#E81E63", "#3F51B5")) +
          labs(y = "Normalized expression", title = gene) + 
          theme(plot.title = element_text(hjust = 0.5, size = 18, face = "italic")) + 
          stat.test)
  dev.off()
}
