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



# FUNCTION: GETTING LOQ ####
#this function gets the limit of quantitation (a metric defined by NanoString) which is the geometric mean * standard devation^2 of the negative probes (each ROI has a negative probe control). 

LOQ <- function(a) {
  negprobeindex <- which(rownames(a) %in% "NegProbe-WTX")
  negprobes <- a[negprobeindex, ]
  gmean <- geomean(as.numeric(negprobes))
  stdev <- geosd(as.numeric(negprobes))
  loq <- gmean*stdev^2
  return(loq)
}

# FUNCTION: APPLING LOQ FILTER TO DATA ####
# This function filters the genes that are above the loq in a user defined number of ROIs (min_rois)

loqfilter <- function(data, min_rois) {
  templist = list()
  negprobeindex <- which(rownames(data) %in% "NegProbe-WTX")
  negprobes <- data[negprobeindex, ]
  gmean <- geomean(as.numeric(negprobes))
  stdev <- geosd(as.numeric(negprobes))
  loq <- gmean*stdev^2
  for(x in colnames(data)){
    temp <- data %>% dplyr::filter(eval(parse(text=x)) > loq)
    templist <- append(rownames(temp), templist)
  }
  tempunlist <- unlist(templist)
  table1 <- table(tempunlist)
  namelist <- names(which(table1 > min_rois))
  index1 <- which(rownames(data) %in% namelist)
  result <- data[index1,]
  return(result)
}

# FUNCTION: volcanoPlot() ####
# This function generates volcano plots and from the output of differential expression analysis.
# input = limma/TopTable output, title_1/2 = string for the groups in order that they appear in the contrast (i.e. lvi vs stroma), batch = string
volcanoPlot <- function(input, title_1, title_2, batch, upcolor = "#619CFF", downcolor = "#F8766D", directory = "/outputs/") {
  # 
  # Making uppers/mids/downers
  data <- input %>%
    mutate(neg.log.padj = -1*log10(adj.P.Val)) %>% 
    filter(!is.na(neg.log.padj)) %>%
    mutate(topstat = neg.log.padj * abs(logFC * 2))
  #data$rn <- rownames(data)
  #data
  Uppers <- data %>% 
    filter((logFC >= 0.58) & (neg.log.padj >= 1.3)) %>% 
    mutate(Group = 'Uppers')
  Mids <- data %>% 
    filter(((logFC < 0.58) & (logFC >-0.58)) | (neg.log.padj < 1.3)) %>%
    mutate(Group = 'Mids')
  Lowers <- data %>% 
    filter((logFC <= -0.58) & (neg.log.padj >= 1.3)) %>%
    mutate(Group = 'Lowers')
  
  # Topstat for labeling volcano plot
  inUpperstopstat = character()
  inLowerstopstat = character()
  inMidstopstat = character ()
  
  
  #Uppers
  if(nrow(Uppers) >= 1){
    for(x in 1:nrow(Uppers)) {
      if(round(Uppers$topstat[x]) %in% round(sort(Uppers$topstat, decreasing = TRUE)[1:6])) {
        inUpperstopstat <- append(inUpperstopstat, TRUE) }
      else{
        inUpperstopstat <- append(inUpperstopstat, FALSE)
      }
    }
  } else {
    print(Uppers)
  }
  
  #Lowers
  for(x in 1:nrow(Mids)){
    inMidstopstat <- append(inMidstopstat, FALSE)
  }
  #Lowers
  if(nrow(Lowers) >= 1) {
    for(x in 1:nrow(Lowers)) {
      if(round(Lowers$topstat[x]) %in% round(sort(Lowers$topstat, decreasing = TRUE)[1:6])) {
        inLowerstopstat <- append(inLowerstopstat, TRUE) }
      else{
        inLowerstopstat <- append(inLowerstopstat, FALSE)
      }
    }
  }
  else{
    print(Lowers)
  }
  
  print(Mids)
  
  midupp <- append(inUpperstopstat, inMidstopstat)
  tops <- append(midupp, inLowerstopstat)
  print(length(tops))
  Volcano_groups <- rbind(Uppers,Mids)
  Volcano_groups2 <- rbind(Volcano_groups,Lowers)
  Volcano_groups3 <- cbind(Volcano_groups2, tops)
  setDT(Volcano_groups3, keep.rownames = TRUE)[]
  #Making DEGs list
  
  degs <- rbind(Uppers, Lowers)
  degs <- degs %>% mutate(statistic = (logFC/abs(logFC))*neg.log.padj)
  
  degs_ordered <- degs[order(degs$statistic, decreasing = TRUE),]
  setDT(degs_ordered, keep.rownames = TRUE)[]
  degs_ordered %>% dplyr::select(c('rn', 'statistic')) -> degsforgsea
  degs_ordered %>% dplyr::select(c('rn', 'logFC', 'AveExpr', "P.Value", "adj.P.Val", "neg.log.padj", "Group")) -> degslist
  colnames(degslist)[1] <- "Genes"
  
  # making gsea lists for all genes ranked either by LFC or "statistic
  allgenes <- rbind(Uppers, Mids, Lowers)
  allgenes <- allgenes %>% mutate(statistic = (logFC/abs(logFC))*neg.log.padj)
  setDT(allgenes, keep.rownames = TRUE)[]
  allgenes_ranked_by_statistic <- allgenes[order(allgenes$statistic, decreasing = TRUE),]
  allgenes_ranked_by_LFC <- allgenes[order(allgenes$logFC, decreasing = TRUE),]
  allgenes_statistic_gsea <- allgenes_ranked_by_statistic %>% dplyr::select(c('rn', 'statistic'))
  allgenes_LFC_gsea <- allgenes_ranked_by_LFC %>% dplyr::select(c('rn', 'logFC'))
  
  
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
  
  volcano_plot  <- ggplot(Volcano_groups3, mapping = aes(x = logFC, y=neg.log.padj, label = rn)) +
    geom_point(mapping = aes(color = Group), size = 1) +
    xlab(expression(paste(log[2], " fold change"))) +
    ylab(expression(paste(-log[10],"(q)"))) +  
    geom_text_repel(data = . %>% 
                      mutate(label = ifelse(Group %in% c("Uppers", "Lowers") & tops == TRUE,
                                            rn, "")),
                    aes(label = label), 
                    show.legend = FALSE,
                    max.overlaps = Inf,
                    force = 10,                    # Increased from 2 to 10
                    force_pull = 2,               # Added: pulls labels toward points
                    box.padding = 0.5,            # Increased from 0.2 to 0.5
                    point.padding = 0.3,          # Added: padding around points
                    min.segment.length = 0.1,     # Added: minimum segment length
                    seed = 42,                    # Added: for reproducible positioning
                    fontface = "italic") +
    ggtitle(paste0(title_2," versus ", title_1)) +
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
    geom_vline(xintercept = c(-0.58,0.58))
  
  volcano_plot
  ggsave(file = paste0(directory, batch,"_",title_1,"vs",title_2,"_Batch_Corrected.pdf"), device = cairo_pdf)
}

# FUNCTION: genes_start_with() ####
# this function is executed within makeplotfolder()
genes_start_with <- function(data, starting_string) {
  genes <- as.data.frame(t(data)) %>% select(starts_with(paste0(starting_string))) %>% colnames()
  return(genes)
}


# FUNCTION: makeplotfolder() ####
# this function makes folder that starts with the prefix of genes (i.e. "MUC"), and then fills the folder with barplots of genes starting with "MUC" (i.e. "MUC13", "MUCL3") 
makeplotfolder <- function(data, starting_string, makedirectory){
  yee <- genes_start_with(data, starting_string)
  for(x in yee){
    gene_column_plot(data = data, gene = x, mkdirectory = makedirectory, directoryname = starting_string)
  }
}



