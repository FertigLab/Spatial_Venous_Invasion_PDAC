library(tidyverse)
library(org.Hs.eg.db)
library(purrr)
library(biomaRt)
library(strex)
library(DESeq2)
library(projectR)
library(limma)
library(edgeR)
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

# reading in custom functions
source("scripts/00_Custom_Functions_CTC.R")

# reading in the raw counts from the CTC bulk RNA seq data
rawcounts <- read.csv("inputs/GSE144561_rawCountsAllsamples.txt", sep = "\t")

# renaming the samples
samplenames <- c("HD1",	
                 "HD2",
                 "HD3",
                 "HD4",	
                 "HD5",
                 "HD6",	
                 "HD7",	
                 "HD8",	
                 "HD9",	
                 "HD10",	
                 "HD11",	
                 "HD12",	
                 "HD13",	
                 "HD14",	
                 "HD15",	
                 "HD16",	
                 "HD17",	
                 "HD18",	
                 "HD19",	
                 "HD20",	
                 "HD21",	
                 "locPDAC01",	
                 "locPDAC02",
                 "locPDAC03",	
                 "locPDAC04",	
                 "locPDAC05",	
                 "locPDAC06",	
                 "locPDAC07",
                 "locPDAC08",
                 "locPDAC09",	
                 "locPDAC10",	
                 "locPDAC11",	
                 "locPDAC12",	
                 "locPDAC13",	
                 "locPDAC14",	
                 "locPDAC15",	
                 "locPDAC16",	
                 "locPDAC17",	
                 "metPDAC01",	
                 "metPDAC02",	
                 "metPDAC03",	
                 "metPDAC04",
                 "metPDAC05",
                 "metPDAC06",	
                 "metPDAC07",	
                 "metPDAC08",	
                 "metPDAC09",	
                 "metPDAC10",	
                 "metPDAC11",	
                 "metPDAC12",	
                 "metPDAC13",	
                 "metPDAC14",	
                 "metPDAC15",	
                 "metPDAC16",	
                 "metPDAC17",	
                 "metPDAC18",	
                 "SU2C01",	
                 "SU2C02",	
                 "SU2C03",	
                 "SU2C04",	
                 "SU2C05",	
                 "SU2C06",	
                 "SU2C07",	
                 "SU2C08",	
                 "SU2C09",	
                 "SU2C10",	
                 "SU2C11",	
                 "SU2C12",	
                 "SU2C13",	
                 "SU2C14",	
                 "SU2C15",	
                 "SU2C16",	
                 "SU2C17",	
                 "SU2C18",	
                 "SU2C19",	
                 "SU2C20",	
                 "SU2C21",	
                 "SU2C22",	
                 "SU2C23",	
                 "SU2C24",	
                 "SU2C25")

# making a phenotable
sampletype <- gsub('[[:digit:]]+', '', samplenames)
sampletypecombined <- gsub('SUC', 'LOCAL_DISEASE', sampletype)
sampletypecombined <- gsub('locPDAC', 'LOCAL_DISEASE', sampletypecombined)
sampletypecombined <- gsub('metPDAC', 'METASTATIC_DISEASE', sampletypecombined)
sampletypecombined <- gsub('HD', 'NORMAL_DONOR', sampletypecombined)

# reading in the cogaps object
cogaps <- readRDS("inputs/batch_merged_cogapsP4.rds")

# subsetting the CTC data for the genes present in the processed geomx data
geomxgenes <- rownames(cogaps@loadingStdDev)
newrawcounts <- na.omit(rawcounts[geomxgenes,])
tingpheno <- as.data.frame(cbind(samplenames, sampletype, sampletypecombined))
rownames(tingpheno) <- colnames(rawcounts)

## DESeq preprocessing with vst
dds <- DESeqDataSetFromMatrix(countData = newrawcounts,
                              colData = tingpheno,
                              design= ~0 +sampletype)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

# doing DESeq2 differential expression analysis
resloc <- results(dds, contrast = list("sampletypelocPDAC", "sampletypeHD"))
resloc <- as.data.frame(resloc)

# getting the significantly upregulated and downregulated genes between circulating tumor cells and healthy donors
sigsloc <- resloc %>% filter(padj < 0.05 & (log2FoldChange >= 0.58))
downloc <- resloc %>% filter(padj < 0.05 & (log2FoldChange <= -0.58))

# FIGURE S11A ---- volcano plot for circulating tumor cells vs healthy donors ####
volcanoPlot(resloc, "CTC", "HD", "deseq2", upcolor = "#FF5722", downcolor = "#000000", directory = "outputs/plots_tables_objects/")
#####

# reading in the pattern markers
jp4thresh <- readRDS('outputs/thresholded_patterns.rds')
jp4thresh$PatternMarkers[[1]] -> pat1markers
jp4thresh$PatternMarkers[[2]] -> pat2markers
jp4thresh$PatternMarkers[[3]] -> pat3markers
jp4thresh$PatternMarkers[[4]] -> pat4markers

pat1markers <- intersect(pat1markers, rownames(resloc))
pat3markers <- intersect(pat3markers, rownames(resloc))


# getting the pattern 1 markers which are significantly upregulated in CTC relative to healthy donor
up1 <- intersect(pat1markers, rownames(sigsloc))

# getting the pattern 3 markers which are significantly upregulated in CTC relative to healthy donors
up3 <- intersect(pat3markers, rownames(sigsloc))

# getting the pattern 1 markers which are significantly downregulated in CTC relative to healthy donors
down1 <- intersect(pat1markers, rownames(downloc))

# getting the pattern 3 markers which are significantly downregulated in CTC relative to healthy donors
down3 <- intersect(pat3markers, rownames(downloc))


# Performing fisher's exact test for the bar plot

# getting number of genes in each category
pat1nonHD <- length(pat1markers) - length(down1)
pat1nonCTC <- length(pat1markers) - length(up1)
pat1nonsig <- length(pat1markers) - length(down1) - length(up1)

pat3nonHD <- length(pat3markers) - length(down3)
pat3nonCTC <- length(pat3markers) - length(up3)
pat3nonsig <- length(pat3markers) - length(down3) - length(up3)

# making the outcomes table for the odds ratios/fisher's test
HD_outcomes <- c(rep(1, length(down1)), rep(0, pat1nonHD), rep(1, length(down3)), rep(0, pat3nonHD))
CTC_outcomes <- c(rep(1, length(up1)), rep(0, pat1nonCTC), rep(1, length(up3)), rep(0, pat3nonCTC))

pat1_vs_pat3_HDaries <- table(
  Population = c(rep("pat1", length(pat1markers)), rep("pat3", length(pat3markers))),
  Outcome = HD_outcomes)

pat1_vs_pat3_CTCs <- table(
  Population = c(rep("pat1", length(pat1markers)), rep("pat3", length(pat3markers))),
  Outcome = CTC_outcomes)

# the odds ratios and p values for the HDary genes
theodds_HDs <- oddsratio(pat1_vs_pat3_HDaries, method = "fisher")

# the odds ratios and p values for the CTC genes
theodds_CTCs <- oddsratio(pat1_vs_pat3_CTCs, method = "fisher")

the_df <- data.frame(matrix(NA, nrow = 256, ncol = 1))
the_df$labeled_geomx <- c(rep("Pattern_1", length(pat1markers)), rep("Pattern_3", length(VIgenes)))
the_df$labeled_CTC <- c(rep("HDary", length(down1)), rep("neither", pat1nonsig), rep("CTC", length(up1)), rep("HD", length(down3)), rep("neither", pat3nonsig), rep("CTC", length(up3)))
the_df$count <- 1



# getting all of the significantly expressed genes 
allsigs <- resloc %>% filter(padj < 0.05 & (log2FoldChange >= 0.58 | log2FoldChange <= -0.58)) 

# ranking the DEGs by log fold change
the_rank <- allsigs$log2FoldChange
names(the_rank) <- rownames(allsigs)
lfcrank <- sort(the_rank, decreasing = T)

# making data table for comparing the overlaps between the CoGAPS patterns and the CTC vs HD DEGs
allthegeomxgenes <- c(pat3markers,pat1markers)
vi_or_pdac <- c(rep("VI", length(pat3markers)), rep("PDAC", length(pat1markers)))
ctc_or_hd <- rep("neither", length(allthegeomxgenes))
ups <- lfcrank > 0
ups <- ups[which(ups == TRUE)]
downs <- lfcrank < 0
downs <- downs[which(downs == TRUE)]
ctc_or_hd[which(allthegeomxgenes %in% names(ups))] <- "CTCs"
ctc_or_hd[which(allthegeomxgenes %in% names(downs))] <- "Healthy donor"
the_df <- as.data.frame(cbind(vi_or_pdac, allthegeomxgenes, ctc_or_hd))
the_df$count <- 1
the_df$ctc_or_hd <- factor(the_df$ctc_or_hd, levels = c("CTCs", "neither", "Healthy donor"))

# FIGURE 7D ---- pattern 1 and 3 markers which are either elevated in CTC or healthy donor ####
pdf("outputs/plots_tables_objects/barplot_patternmarkers_CTC.pdf")
ggplot(the_df) + 
  geom_bar(aes(fill=ctc_or_hd, y=count, x=vi_or_pdac), position="stack", stat="identity") + 
  scale_fill_manual(name = "CTC RNA-Seq", labels = c("Significantly upregulated in CTC", "Not signficantly upregulated", "Significantly upregulated in HD"), values = c("#FF5722", "#9C9C9C", "#000000")) +
  scale_x_discrete(name = "GeoMx", labels = c("Pattern 1 Markers", "Pattern 3 Markers")) + 
  ylab("Number of Significantly Upregulated Genes")
dev.off()
#####

# Making the gene box plots

# reading in the normalized data from the CTC paper
prenormalized <- read.csv("inputs/GSE144561_DESeq2_NormalizedCountsForAllSamples.txt", sep = "\t")

# log transforming the data
prenormalized <- log10(prenormalized + 1)

# filtering for the HD and localized + treatment naive PDAC samples
prenormalized <- prenormalized[,1:38]

# scaling the data
datScale <- t(apply(prenormalized,1,scale))

# transposing the data and combining it with the metadata
qcDat <- data.frame(tingpheno$sampletype[1:38],t(datScale))

# changing the colnames and groupnames
colnames(qcDat)[1] <- "Group"
qcDat$Group[which(qcDat$Group == "locPDAC")] <- "CTC"

# FIGURE 7E, S11B-S11C ####
# making box plots for the pattern 3 markers which are significantly upregulated in CTC relative to healthy donors. Adding stats too (padj from DESEq2)
ctc_hd <- data.frame("Genes" = rownames(sigsloc), "group_1" = "CTC", "group_2" = "HD", "p" = sigsloc$padj)
mygenes <- intersect(pat3markers, rownames(sigsloc))
for(x in mygenes){
  data_for_p_val_manual <- ctc_hd %>% filter(Genes == x) %>% dplyr::select("group_1", "group_2", "p")
  colnames(data_for_p_val_manual) <- c("group1", "group2", "p.adj")
  datman <- tibble(data_for_p_val_manual)
  
  # Convert p-values to asterisks
  datman$p.signif <- sapply(datman$p.adj, p_to_asterisk)
  
  yposition <- max(qcDat[,x]) * 1.05
  stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
  pdf(paste0("outputs/geneplots/ctc_pattern_3/",x,".pdf"), width = 4, height = 6)
  print(ggplot(qcDat, mapping = aes_string(x = "Group", y=x, color = "Group")) + geom_boxplot() +
          geom_jitter(shape=16, position=position_jitter(0.2))+
          scale_color_manual(name = "", labels = c("CTC", "HD"), values = c("#FF5722", "#000000")) +
          labs(y = "Normalized expression", title = x) + 
          theme(plot.title = element_text(hjust = 0.5, size = 18, face = "italic")) + 
          stat.test)
  dev.off()
}

# making box plots for the pattern 1 markers which are significantly upregulated in CTC relative to healthy donors. Adding stats too (padj from DESEq2)
mygenes <- intersect(pat1markers, rownames(sigsloc))
for(x in mygenes){
  data_for_p_val_manual <- ctc_hd %>% filter(Genes == x) %>% dplyr::select("group_1", "group_2", "p")
  colnames(data_for_p_val_manual) <- c("group1", "group2", "p.adj")
  datman <- tibble(data_for_p_val_manual)
  
  # Convert p-values to asterisks
  datman$p.signif <- sapply(datman$p.adj, p_to_asterisk)
  
  yposition <- max(qcDat[,x]) * 1.05
  stat.test <- stat_pvalue_manual(datman, y.position = yposition, step.increase = 0.1, label = "p.signif", hide.ns = F, size = 6)
  pdf(paste0("outputs/geneplots/ctc_pattern_1/",x,".pdf"), width = 4, height = 6)
  print(ggplot(qcDat, mapping = aes_string(x = "Group", y=x, color = "Group")) + geom_boxplot() +
          geom_jitter(shape=16, position=position_jitter(0.2))+ 
          scale_color_manual(name = "", labels = c("CTC", "HD"), values = c("#FF5722", "#000000")) +
          labs(y = "Normalized expression", title = x) + 
          theme(plot.title = element_text(hjust = 0.5, size = 18, face = "italic")) + 
          stat.test)
  dev.off()
}

#####