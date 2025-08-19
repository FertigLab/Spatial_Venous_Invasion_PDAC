

library(ComplexHeatmap)
library(knitr)
#library(circlize)
#library(colorspace)
#library(GetoptLong)
#library(psych)
#library(FSA)
library(tidyverse)
library(sva)
library(limma)
library(Biobase)
library(patchwork)
library(data.table)
library(ggrepel)
library(DESeq2)
library(DT)
library(EnhancedVolcano)
#library(dplyr)
#library(ggforce)
library(ggplot2)
library(webshot)

# FUNCTION: TURNING P VALUES INTO ASTERISKS
p_to_asterisk <- function(p_val) {
  if (p_val < 0.00001) {
    return("*****")
  } else if (p_val < 0.0001) {
    return("****")
  } else if (p_val < 0.001) {
    return("***")
  } else if (p_val < 0.01) {
    return("**")
  } else if (p_val < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}




# FUNCTION: volcanoPlot() ####
# This function generates volcano plots and from the output of differential expression analysis FOR DESEQ2 NOT LIMMA
# input = limma/TopTable output, title_1/2 = string for the groups in order that they appear in the contrast (i.e. lvi vs stroma), batch = string
volcanoPlot <- function(input, title_1, title_2, batch, upcolor = "#619CFF", downcolor = "#F8766D",  lfcthresh = 0.58, directory = "/Users/alexanderbell/OneDrive - Johns Hopkins/GeoMx_Chianchiano_Project/March_2023_Ting_CTC/outputs/") {
  # Making uppers/mids/downers
  data <- input %>%
    mutate(neg.log.padj = -1*log10(padj)) %>% 
    filter(!is.na(neg.log.padj)) %>%
    mutate(topstat = neg.log.padj * abs(log2FoldChange * 2))
  #data$rn <- rownames(data)
  #data
  Uppers <- data %>% 
    filter((log2FoldChange >= lfcthresh) & (neg.log.padj >= 1.3)) %>% 
    mutate(Group = 'Uppers')
  Mids <- data %>% 
    filter(((log2FoldChange < lfcthresh) & (log2FoldChange >-lfcthresh)) | (neg.log.padj < 1.3)) %>%
    mutate(Group = 'Mids')
  Lowers <- data %>% 
    filter((log2FoldChange <= -lfcthresh) & (neg.log.padj >= 1.3)) %>%
    mutate(Group = 'Lowers')
  
  # Topstat for labeling volcano plot
  inUpperstopstat = character()
  inLowerstopstat = character()
  inMidstopstat = character ()
  
  
  #Uppers
  for(x in 1:nrow(Uppers)) {
    if(round(Uppers$topstat[x]) %in% round(sort(Uppers$topstat, decreasing = TRUE)[1:6])) {
      inUpperstopstat <- append(inUpperstopstat, TRUE) }
    else{
      inUpperstopstat <- append(inUpperstopstat, FALSE)
    }
  }
  
  #Lowers
  for(x in 1:nrow(Mids)){
    inMidstopstat <- append(inMidstopstat, FALSE)
  }
  for(x in 1:nrow(Lowers)) {
    if(round(Lowers$topstat[x]) %in% round(sort(Lowers$topstat, decreasing = TRUE)[1:6])) {
      inLowerstopstat <- append(inLowerstopstat, TRUE) }
    else{
      inLowerstopstat <- append(inLowerstopstat, FALSE)
    }
  }
  
  midupp <- append(inUpperstopstat, inMidstopstat)
  tops <- append(midupp, inLowerstopstat)
  Volcano_groups <- rbind(Uppers,Mids)
  Volcano_groups2 <- rbind(Volcano_groups,Lowers)
  Volcano_groups3 <- cbind(Volcano_groups2, tops)
  print(Volcano_groups3)
  setDT(Volcano_groups3, keep.rownames = TRUE)[]
  
  #Making DEGs list
  
  degs <- rbind(Uppers, Lowers)
  degs <- degs %>% mutate(statistic = (log2FoldChange/abs(log2FoldChange))*neg.log.padj)
  
  degs_ordered <- degs[order(degs$statistic, decreasing = TRUE),]
  setDT(degs_ordered, keep.rownames = TRUE)[]
  degs_ordered %>% dplyr::select(c('rn', 'statistic')) -> degsforgsea
  degs_ordered %>% dplyr::select(c('rn', 'log2FoldChange', 'baseMean', "pvalue", "padj", "neg.log.padj", "Group")) -> degslist
  colnames(degslist)[1] <- "Genes"
  
  # making gsea lists for all genes ranked either by LFC or "statistic
  allgenes <- rbind(Uppers, Mids, Lowers)
  allgenes <- allgenes %>% mutate(statistic = (log2FoldChange/abs(log2FoldChange))*neg.log.padj)
  setDT(allgenes, keep.rownames = TRUE)[]
  allgenes_ranked_by_statistic <- allgenes[order(allgenes$statistic, decreasing = TRUE),]
  allgenes_ranked_by_LFC <- allgenes[order(allgenes$log2FoldChange, decreasing = TRUE),]
  allgenes_statistic_gsea <- allgenes_ranked_by_statistic %>% dplyr::select(c('rn', 'statistic'))
  allgenes_LFC_gsea <- allgenes_ranked_by_LFC %>% dplyr::select(c('rn', 'statistic'))
  
  
  # all genes
  write.table(Volcano_groups3, file = paste0(directory, batch,"_", title_1, "vs", title_2, "_all_genes.csv"), row.names = TRUE, col.names = TRUE, sep = ',')
  
  # just degs, gsea formatted
  write.table(degsforgsea, file = paste0(directory, batch,"_", title_1, "vs", title_2, "_gseadegs.csv"), row.names = TRUE, col.names = TRUE, sep = ',')
  
  # just degs
  write.table(degslist, file = paste0(directory, batch, "_", title_1, "vs", title_2, "_DEGS.csv"), row.names = TRUE, col.names = TRUE, sep = ',')
  
  # all genes, gsea formatted, ranked by "statistic" (Aka the neg.log.fold change and direction)
  write.table(allgenes_statistic_gsea, file = paste0(directory, batch, "_", title_1, "vs", title_2, "_gsea_allgenes_rankedstatistic.csv"), row.names = TRUE, col.names = TRUE, sep = ',')
  
  # all genes, gsea formatted, ranked by LFC
  write.table(allgenes_LFC_gsea, file = paste0(directory, batch, "_", title_1, "vs", title_2, "_gsea_allgenes_rankedLFC.csv"), row.names = TRUE, col.names = TRUE, sep = ',')
  
  volcano_plot  <- ggplot(Volcano_groups3, mapping = aes(x = log2FoldChange, y=neg.log.padj, label = rn)) +
    geom_point(mapping = aes(color = Group), size = 1) +
    xlab(expression(paste(log[2], " fold change"))) +
    ylab(expression(paste(-log[10],"(q)"))) +  
    geom_text_repel(data = . %>% 
                      mutate(label = ifelse(Group %in% c("Uppers", "Lowers") & tops == TRUE,
                                            rn, "")),
                    aes(label = label), 
                    show.legend = FALSE,
                    max.overlaps = Inf,
                    force = 2,
                    box.padding = 0.2,
                    fontface = "italic") +
    ggtitle(paste0(title_1," versus", title_2)) +
    scale_color_manual(name = "Groups", labels = c(paste0("Increased in ", title_2, ": ", nrow(Lowers)), paste0("Unchanged: ", nrow(Mids)), paste0("Increased in ", title_1, ": ", nrow(Uppers))), 
                       values = c(downcolor, "#CBCBCB" , upcolor )) +
    theme_bw()+
    theme(text = element_text(family = "Arial", color = "#2C2C2E"),
          plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
          legend.title = element_blank(),
          axis.text=element_text(size=10), 
          axis.title=element_text(size=10,face="bold"),
          legend.text = element_text(size=10),
          legend.position=c(1.02,0.5),
          legend.justification=c(0, 1), 
          legend.key.width=unit(1, "lines"), 
          legend.key.height=unit(1, "lines"), 
          plot.margin = unit(c(1, 13, 0.5, 0.5), "lines"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    geom_hline(yintercept = 1.3) +
    geom_vline(xintercept = c(-lfcthresh,lfcthresh))
  
  volcano_plot
  ggsave(file = paste0(directory,batch,"_",title_1,"vs",title_2,".pdf"), device = cairo_pdf)
}