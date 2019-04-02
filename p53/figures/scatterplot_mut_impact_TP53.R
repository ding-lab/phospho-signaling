# Yige Wu @ WashU 2019 Jan
## 

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/expression_matrices/expression_shared.R')

library(ggrepel)
library(readxl)
library(ggpubr)

# set variables -----------------------------------------------------------
sig_thres <- 0.05
networkin_thres <- 5
cancers2process <- cancers_sort
amp_thres <- 0.1
del_thres <- -0.1
num_genoalt_thres <- 4
my_comparisons <- list( geneA_mutation = c("geneA_mutation", "control"))
TP53_alteration_sort <- c("TAD1", "TAD2", "proline rich\ndomain", "DBD", "NLS", "tetramerization\ndomain", "regulatory\ndomain", "other", 
                          paste0("TP53", "\n", "deletion"), paste0("TP53", "\n", "amplification"), "control")
position_range <- c(170, 300)
hotspot <- c(175, 248, 273)
hotspot_proximal <- as.vector(sapply(hotspot, function(i) i+(-5):5))

# inputs ------------------------------------------------------------------


# decide what to plot -----------------------------------------------------
# pairs2plot <- paste0("TP53", ":", "ESR1", ":", "PRO")
pairs2plot <- paste0("TP53", ":", "ESR1", ":", "RNA")
# pairs2plot <- paste0("TP53", ":", "TP53", ":", "RNA")
# pairs2plot <- paste0("TP53", ":", "TP53", ":", "PRO")
# pairs2plot <- paste0("TP53", ":", "CDKN1A", ":", "RNA")
# pairs2plot <- paste0("TP53", ":", "CDKN1A", ":", "PRO")
# pairs2plot <- paste0("TP53", ":", "MDM2", ":", "RNA")
# pairs2plot <- paste0("TP53", ":", "CDKN1A", ":", "PRO")

pairs2plot

# plot per pair all cancers-----------------------------------------------------------
for (pair in pairs2plot) {
  geneA <- str_split(string = pair, pattern = ":")[[1]][1]
  geneB <- str_split(string = pair, pattern = ":")[[1]][2]
  phosphosite <- str_split(string = pair, pattern = ":")[[1]][3]
  fn = paste0(makeOutDir(resultD = resultD), geneA, "_", geneB, "_", phosphosite, "_by_domain_linear.pdf")
  
  # for (cancer in c("BRCA")) {
  if (!file.exists(fn)) {
    sup_tab <- NULL
    
    for (cancer in c("BRCA", "CO", "UCEC", "OV", "CCRCC", "LIHC")) {
      
      # input data first because different for each cancer type --------------------------------------------------------------
      ## input maf
      maf <- loadMaf(cancer = cancer, maf_files = maf_files)
      maf <- maf[maf$Hugo_Symbol == geneA,]
      
      ## input CNA matrix
      if (geneA == "TP53") {
        cna_tab <- loadTP53Deletion(cancer = cancer)
      } else {
        cna_tab <- loadCNAstatus(cancer = cancer)
      }
      cna_tab <- cna_tab[cna_tab$gene %in% c(geneA), ]
      
      ## load RNA
      rna_tab <- loadRNA(cancer = cancer)
      rna_tab <- rna_tab[rna_tab$gene %in% c(geneB),]
      
      ## input protein data
      if (cancer %in% c("BRCA", "OV", "CO")) {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
        
      } else if (cancer == "UCEC") {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
        
      } else if (cancer == "CCRCC") {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
      } else if (cancer == "LIHC") {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
        
      }
      pro_tab <- pro_tab[pro_tab$Gene %in% c(geneB),]
      pho_tab <- pho_tab[pho_tab$Gene == geneB & pho_tab$Phosphosite == phosphosite,]
      
      # make the annotation columns for each sample -----------------------------
      partIDs <- colnames(pro_tab)[!(colnames(pro_tab) %in% c("Gene", "Phosphosite", "Peptide_ID"))]
      col_anno <- data.frame(partID = partIDs)
      
      if (nrow(maf) > 0){
        maf$partID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
        maf.m <- melt(maf[, c("partID", "HGVSp_Short")], id.vars = c("partID"))
        maf.m %>% head()
        maf.m <- maf.m[, c("partID", "value")]
        colnames(maf.m) <- c("partID", "HGVSp_Short")
        
        col_anno <- merge(col_anno, maf.m, by = c("partID"), all.x = T)
      } else {
        print("no mutation!")
      }
      
      ## CNA needs to show both geneA and geneB; levels: amplification, deletion, neutral
      if (nrow(cna_tab) > 0) {
        cna_tab.m <- melt(cna_tab, id.vars = "gene")
        colnames(cna_tab.m) <- c("gene", "partID", "CNA")
        
        for (gene in unique(cna_tab.m$gene[cna_tab.m$CNA != "neutral"])) {
          cna_mat2merge <- cna_tab.m[cna_tab.m$gene == gene, c("partID", "CNA")]
          colnames(cna_mat2merge) <- c("partID", paste0("CNA.", gene))
          col_anno <- merge(col_anno, cna_mat2merge, by = c("partID"), all.x = T)
        }
      } else {
        print("no CNA!")
      }
      
      ## order samples
      col_anno %>% head()
      
      for (gene in c(geneA, geneB)) {
        if (paste0("CNA.", gene) %in% colnames(col_anno)) {
          col_anno <- col_anno[order(col_anno[, paste0("CNA.", gene)], decreasing = T),]
          ann_colors[[paste0("CNA.", gene)]] <-  c(amplification = "#E41A1C", deletion = "#377EB8", "neutral" = "grey")
        }
        if (paste0("mutation.", gene) %in% colnames(col_anno)) {
          col_anno <- col_anno[order(col_anno[, paste0("mutation.", gene)], decreasing = T),]
          ann_colors[[paste0("mutation.", gene)]] <- c(missense = "#E41A1C", truncation = "#377EB8", wild_type = "white", "missense&truncation" = "#6A3D9A", other_mutation = "#FF7F00", silent = "#33A02C")
        }
      }
      
      # make the matrix of values showing in heatmap ----------------------------
      sup_tab_can <- NULL
      if (nrow(pro_tab) > 0) {
        pro_tab.m <- melt(pro_tab, id.vars = "Gene")
        pro_tab.m %>% head()
        colnames(pro_tab.m) <- c("Gene", "partID", "exp_value")
        pro_tab.m$Phosphosite <- "PRO"
        sup_tab_can <- rbind(sup_tab_can, pro_tab.m[,c("Phosphosite", "partID", "exp_value")])
      }
      
      rna_tab.m <- melt(rna_tab, id.vars = "gene")
      colnames(rna_tab.m) <- c("Gene", "partID", "exp_value")
      rna_tab.m$Phosphosite <- "RNA"
      if (nrow(rna_tab.m) > 0) {
        sup_tab_can <- rbind(sup_tab_can, rna_tab.m[,c("Phosphosite", "partID", "exp_value")])
      }
      # sup_tab_can$id_row <- paste0(sup_tab_can$Gene, "_", sup_tab_can$Phosphosite)
      # sup_tab <- sup_tab[sup_tab$id_row == paste0(geneA, "_", "PRO"),]
      
      if (!is.null(sup_tab_can)) {
        sup_tab_can$exp_value <- as.numeric(as.vector(sup_tab_can$exp_value))
        sup_tab_can <- unique(sup_tab_can)
        
        # input sample classification ---------------------------------------------
        sample_class <- fread(input = paste0("./Ding_Lab/Projects_Current/TP53_shared_data/resources/TP53_samples_classification/", cancer, "_tp53_mut_classification_table.txt"), data.table = F)
        if (cancer == "LIHC") {
          meta_tab <- read_excel("./Ding_Lab/Projects_Current/TP53_shared_data/resources/TP53_samples_classification/HCC Clinical information and follow-up information.xlsx", 
                                 sheet = "Table S1 Clinical and follow-up")
          meta_tab <- data.frame(meta_tab)
          meta_tab$partID <- paste0("LIHC", meta_tab$Case.ID)
          rownames(meta_tab) <- as.vector(meta_tab$partID)
          meta_tab$sampID.tumor <- str_split_fixed(string = meta_tab$Tumor..T..sample.ID, pattern = "T", 2)[,2]
          meta_tab$sampID.normal <- str_split_fixed(string = meta_tab$Adjacent.Non.tumor.liver.tissue..N..sample.ID, pattern = "N", 2)[,2]
          sample_class$partID <- sapply(X = sample_class$V1, FUN = function(sampID, meta_tab) {
            partID <- unlist(meta_tab$partID[meta_tab$sampID.normal == sampID | meta_tab$sampID.tumor == sampID])[1]
            return(partID)
          }, meta_tab = meta_tab)
        } else {
          colnames(sample_class)[1] <- "partID"
          
        }
        col_anno <- merge(col_anno, sample_class[, c("partID", "Classification_complex")], all.x = T)
        sup_tab_can <- merge(sup_tab_can, col_anno, by = c("partID"), all.x = T)
        sup_tab_can$cancer <- cancer
        sup_tab <- rbind(sup_tab, sup_tab_can)
      }
      
    }
  }
  
  # format the data frame to plot -------------------------------------------
  sup_tab2plot <- sup_tab[sup_tab$Phosphosite == phosphosite,]
  tab2p <- sup_tab2plot
  tab2p$position <- str_split_fixed(string = str_split_fixed(string = tab2p$HGVSp_Short, pattern = 'p.[A-Z]', 2)[,2], pattern = '\\*|[A-Z]', n = 2)[,1]
  tab2p$position <- str_split_fixed(string = tab2p$position, pattern = '\\_', n = 2)[,1]
  tab2p$position <- str_split_fixed(string = tab2p$position, pattern = '[a-z]', n = 2)[,1]
  tab2p$position <- as.numeric(as.vector(tab2p$position))
  tab2p$position_y <- tab2p$position - position_range[1]
  tab2p$position_y[grepl(x = tab2p$Classification_complex, pattern = "Truncation")] <- (-20)
  tab2p$position_y[tab2p$CNA.TP53 == "deletion" & is.na(tab2p$position)] <- (-40)
  tab2p$position_y[tab2p$CNA.TP53 == "neutral" & is.na(tab2p$position)] <- (-60)
  
  tab2p$x <- as.vector(tab2p$exp_value)
  tab2p$y <- tab2p$position_y
  tab2p <- tab2p[!is.na(tab2p$cancer),]
  which(is.na(tab2p$y))
  tab2p <- tab2p[(tab2p$position >= position_range[1] & tab2p$position <= position_range[2] & (tab2p$position > position_range[1])) | is.na(tab2p$position) ,]
  tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
  
  tab2p
  # add control average values ----------------------------------------------
  control_exp_med <- sapply(unique(tab2p$cancer), FUN = function(cancer, tab2p) {
    tab2p_can <- tab2p[tab2p$cancer == cancer & tab2p$Classification_complex == "WT",]
    control_exp_med <- median(tab2p_can$exp_value, na.rm = T)
    return(control_exp_med)
  }, tab2p = tab2p)
  control_exp_meds <- data.frame(control_exp_med = control_exp_med, cancer = names(control_exp_med))
  tab2p <- merge(tab2p, control_exp_meds, by = c("cancer"), all.x = T)
  ## get the high quantile of expression values
  exp_high <- sapply(unique(tab2p$cancer), FUN = function(cancer, tab2p) {
    tab2p_can <- tab2p[tab2p$cancer == cancer,]
    exp_high <- quantile(x = tab2p_can$exp_value, probs = 0.9, na.rm = T)
    return(exp_high)
  }, tab2p = tab2p)
  names(exp_high) <-  unique(tab2p$cancer)
  exp_highs <- data.frame(exp_high = exp_high, cancer = unique(tab2p$cancer))
  tab2p <- merge(tab2p, exp_highs, by = c("cancer"), all.x = T)
  ## get the low quantile of expression values
  exp_low <- sapply(unique(tab2p$cancer), FUN = function(cancer, tab2p) {
    tab2p_can <- tab2p[tab2p$cancer == cancer,]
    exp_low <- quantile(x = tab2p_can$exp_value, probs = 0.1, na.rm = T)
    return(exp_low)
  }, tab2p = tab2p)
  names(exp_low) <-  unique(tab2p$cancer)
  exp_lows <- data.frame(exp_low = exp_low, cancer = unique(tab2p$cancer))
  tab2p <- merge(tab2p, exp_lows, by = c("cancer"), all.x = T)
  
  
  # add the position of missense --------------------------------------------
  tab2p$point_fill <- "grey50"
  tab2p$point_fill[tab2p$position_y > 0] <- "orange"
  tab2p$point_fill[(as.vector(tab2p$x) > as.vector(tab2p$exp_high)) & (tab2p$position_y > 0)] <- "red"
  tab2p$point_fill[grepl(x = tab2p$Classification_complex, pattern = "Truncation") & (tab2p$position_y < 0)] <- "blue"
  tab2p$point_fill[tab2p$CNA.TP53 == "deletion" & (tab2p$position_y < 0)] <- "lightblue"
  
  tab2p$cancer <- factor(tab2p$cancer, levels = c("LIHC", "BRCA", "OV", "CO", "UCEC", "CCRCC"))
  
  p = ggplot()
  p = p + geom_point(data = tab2p, mapping = aes(x=x, y=y, fill = point_fill), stroke = 0, alpha = 0.6, size = 2.5, shape = 21)
  p <- p + geom_hline(yintercept = c(hotspot - position_range[1]), linetype = 2, color = "grey50", alpha = 0.5)
  p = p + geom_vline(data = tab2p, mapping = aes(xintercept = control_exp_med), linetype = 2, color = "grey50", alpha = 0.5)
  p = p + geom_vline(data = tab2p, mapping = aes(xintercept = exp_high), linetype = 2, color = set1[1], alpha = 0.5)
  p = p + geom_vline(data = tab2p, mapping = aes(xintercept = exp_low), linetype = 2, color = set1[2], alpha = 0.5)
  p = p + geom_text_repel(data = tab2p[(as.vector(tab2p$x) > as.vector(tab2p$exp_high)) & tab2p$position %in% hotspot_proximal,], 
                          mapping = aes(x=x, y=y, label= as.character(HGVSp_Short)), color = set1[1], segment.color = "black",
                          force = 4, segment.size = 1, segment.alpha = 0.6, size = 4, alpha = 0.8)
  p = p + geom_text_repel(data = tab2p[(as.vector(tab2p$x) < as.vector(tab2p$exp_low)) & tab2p$position %in% hotspot_proximal,], 
                          mapping = aes(x=x, y=y, label= as.character(HGVSp_Short)), color = set1[2], segment.color = "black",
                          force = 4, segment.size = 1, segment.alpha = 0.6, size = 4, alpha = 0.8)
  p <- p +  scale_y_continuous(breaks = c(c(hotspot - position_range[1]), -20, -40, -60), labels = c(as.character(hotspot), "Truncation", "Deletion", "WT"))
  p = p + scale_fill_manual(values = c("red" = set1[1], "blue" = set1[2], "grey50" = "grey50", "orange" = "#FF7F00", "lightblue" = "#A6CEE3"))
  p = p + facet_grid(.~cancer, scales = "free", space = "free_y")
  p = p + theme_nogrid()
  p = p + labs(x = paste0(geneB, "_", phosphosite, " abundance(log2 ratio)"), y = paste0(geneA, " mutation position"))
  p = p + theme(axis.title.x = element_text(size = 15, face = "bold"),
                axis.title.y = element_text(size = 15, face = "bold"),
                axis.text.y = element_text(colour="black", size=8))
  p = p + theme(title = element_text(size = 18, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 8))
  p <- p + theme(panel.spacing.x = unit(0, "lines"),  panel.spacing.y = unit(0, "lines"))
  p <- p + theme(legend.position="none")
  p
  ggsave(file=fn, width = 15, height = 4, useDingbats = F)
  dev.off()
  
}

# plot per pair just BRCA cancers-----------------------------------------------------------
for (pair in pairs2plot) {
  geneA <- str_split(string = pair, pattern = ":")[[1]][1]
  geneB <- str_split(string = pair, pattern = ":")[[1]][2]
  phosphosite <- str_split(string = pair, pattern = ":")[[1]][3]
  
  # for (cancer in c("BRCA")) {
  
  
  for (cancer in c("BRCA")) {
    fn = paste0(makeOutDir(resultD = resultD), cancer, "_", geneA, "_", geneB, "_", phosphosite, "_by_domain_linear.pdf")
    
    if (!file.exists(fn)) {
      sup_tab <- NULL
      # input data first because different for each cancer type --------------------------------------------------------------
      ## input maf
      maf <- loadMaf(cancer = cancer, maf_files = maf_files)
      maf <- maf[maf$Hugo_Symbol == geneA,]
      
      ## input CNA matrix
      cna_tab <- loadCNAstatus(cancer = cancer)
      cna_tab <- cna_tab[cna_tab$gene %in% c(geneA, geneB), ]
      
      ## load RNA
      rna_tab <- loadRNA(cancer = cancer)
      rna_tab <- rna_tab[rna_tab$gene %in% c(geneB),]
      
      ## input protein data
      if (cancer %in% c("BRCA", "OV", "CO")) {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
        
      } else if (cancer == "UCEC") {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
        
      } else if (cancer == "CCRCC") {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
      } else if (cancer == "LIHC") {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
        
      }
      pro_tab <- pro_tab[pro_tab$Gene %in% c(geneB),]
      pho_tab <- pho_tab[pho_tab$Gene == geneB & pho_tab$Phosphosite == phosphosite,]
      
      # make the annotation columns for each sample -----------------------------
      partIDs <- colnames(pro_tab)[!(colnames(pro_tab) %in% c("Gene", "Phosphosite", "Peptide_ID"))]
      col_anno <- data.frame(partID = partIDs)
      
      if (nrow(maf) > 0){
        maf$partID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
        maf.m <- melt(maf[, c("partID", "HGVSp_Short")], id.vars = c("partID"))
        maf.m %>% head()
        maf.m <- maf.m[, c("partID", "value")]
        colnames(maf.m) <- c("partID", "HGVSp_Short")
        
        col_anno <- merge(col_anno, maf.m, by = c("partID"), all.x = T)
      } else {
        print("no mutation!")
      }
      
      ## CNA needs to show both geneA and geneB; levels: amplification, deletion, neutral
      if (nrow(cna_tab) > 0) {
        cna_tab.m <- melt(cna_tab, id.vars = "gene")
        colnames(cna_tab.m) <- c("gene", "partID", "CNA")
        
        for (gene in unique(cna_tab.m$gene[cna_tab.m$CNA != "neutral"])) {
          cna_mat2merge <- cna_tab.m[cna_tab.m$gene == gene, c("partID", "CNA")]
          colnames(cna_mat2merge) <- c("partID", paste0("CNA.", gene))
          col_anno <- merge(col_anno, cna_mat2merge, by = c("partID"), all.x = T)
        }
      } else {
        print("no CNA!")
      }
      
      ## order samples
      col_anno %>% head()
      
      for (gene in c(geneA, geneB)) {
        if (paste0("CNA.", gene) %in% colnames(col_anno)) {
          col_anno <- col_anno[order(col_anno[, paste0("CNA.", gene)], decreasing = T),]
          ann_colors[[paste0("CNA.", gene)]] <-  c(amplification = "#E41A1C", deletion = "#377EB8", "neutral" = "grey")
        }
        if (paste0("mutation.", gene) %in% colnames(col_anno)) {
          col_anno <- col_anno[order(col_anno[, paste0("mutation.", gene)], decreasing = T),]
          ann_colors[[paste0("mutation.", gene)]] <- c(missense = "#E41A1C", truncation = "#377EB8", wild_type = "white", "missense&truncation" = "#6A3D9A", other_mutation = "#FF7F00", silent = "#33A02C")
        }
      }
      
      # make the matrix of values showing in heatmap ----------------------------
      sup_tab_can <- NULL
      if (nrow(pro_tab) > 0) {
        pro_tab.m <- melt(pro_tab, id.vars = "Gene")
        pro_tab.m %>% head()
        colnames(pro_tab.m) <- c("Gene", "partID", "exp_value")
        pro_tab.m$Phosphosite <- "PRO"
        sup_tab_can <- rbind(sup_tab_can, pro_tab.m[,c("Phosphosite", "partID", "exp_value")])
      }
      
      rna_tab.m <- melt(rna_tab, id.vars = "gene")
      colnames(rna_tab.m) <- c("Gene", "partID", "exp_value")
      rna_tab.m$Phosphosite <- "RNA"
      if (nrow(rna_tab.m) > 0) {
        sup_tab_can <- rbind(sup_tab_can, rna_tab.m[,c("Phosphosite", "partID", "exp_value")])
      }
      # sup_tab_can$id_row <- paste0(sup_tab_can$Gene, "_", sup_tab_can$Phosphosite)
      # sup_tab <- sup_tab[sup_tab$id_row == paste0(geneA, "_", "PRO"),]
      
      if (!is.null(sup_tab_can)) {
        sup_tab_can$exp_value <- as.numeric(as.vector(sup_tab_can$exp_value))
        sup_tab_can <- unique(sup_tab_can)
        
        # input sample classification ---------------------------------------------
        sample_class <- fread(input = paste0("./Ding_Lab/Projects_Current/TP53_shared_data/resources/TP53_samples_classification/", cancer, "_tp53_mut_classification_table.txt"), data.table = F)
        if (cancer == "LIHC") {
          meta_tab <- read_excel("./Ding_Lab/Projects_Current/TP53_shared_data/resources/TP53_samples_classification/HCC Clinical information and follow-up information.xlsx", 
                                 sheet = "Table S1 Clinical and follow-up")
          meta_tab <- data.frame(meta_tab)
          meta_tab$partID <- paste0("LIHC", meta_tab$Case.ID)
          rownames(meta_tab) <- as.vector(meta_tab$partID)
          meta_tab$sampID.tumor <- str_split_fixed(string = meta_tab$Tumor..T..sample.ID, pattern = "T", 2)[,2]
          meta_tab$sampID.normal <- str_split_fixed(string = meta_tab$Adjacent.Non.tumor.liver.tissue..N..sample.ID, pattern = "N", 2)[,2]
          sample_class$partID <- sapply(X = sample_class$V1, FUN = function(sampID, meta_tab) {
            partID <- unlist(meta_tab$partID[meta_tab$sampID.normal == sampID | meta_tab$sampID.tumor == sampID])[1]
            return(partID)
          }, meta_tab = meta_tab)
        } else {
          colnames(sample_class)[1] <- "partID"
          
        }
        col_anno <- merge(col_anno, sample_class[, c("partID", "Classification_complex")], all.x = T)
        sup_tab_can <- merge(sup_tab_can, col_anno, by = c("partID"), all.x = T)
        sup_tab_can$cancer <- cancer
        sup_tab <- rbind(sup_tab, sup_tab_can)
      }
      
    }
  }
  
  # format the data frame to plot -------------------------------------------
  sup_tab$subtype <- partID2pam50(patientID_vector = sup_tab$partID, pam50_map = loadPAM50Map())
  sup_tab2plot <- sup_tab[sup_tab$Phosphosite == phosphosite,]
  tab2p <- sup_tab2plot
  tab2p$position <- str_split_fixed(string = str_split_fixed(string = tab2p$HGVSp_Short, pattern = 'p.[A-Z]', 2)[,2], pattern = '\\*|[A-Z]', n = 2)[,1]
  tab2p$position <- str_split_fixed(string = tab2p$position, pattern = '\\_', n = 2)[,1]
  tab2p$position <- str_split_fixed(string = tab2p$position, pattern = '[a-z]', n = 2)[,1]
  tab2p$position <- as.numeric(as.vector(tab2p$position))
  tab2p$position_y <- tab2p$position - position_range[1]
  tab2p$position_y[grepl(x = tab2p$Classification_complex, pattern = "Truncation")] <- (-20)
  tab2p$position_y[tab2p$CNA.TP53 == "deletion" & is.na(tab2p$position)] <- (-40)
  tab2p$position_y[tab2p$CNA.TP53 == "neutral" & is.na(tab2p$position)] <- (-60)
  
  tab2p$x <- as.vector(tab2p$exp_value)
  tab2p$y <- tab2p$position_y
  tab2p <- tab2p[!is.na(tab2p$cancer),]
  which(is.na(tab2p$y))
  tab2p <- tab2p[(tab2p$position >= position_range[1] & tab2p$position <= position_range[2] & (grepl(x = tab2p$Classification_complex, pattern = "Missense"))) | !(grepl(x = tab2p$Classification_complex, pattern = "Missense")) ,]
  tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
  
  tab2p
  # add control average values ----------------------------------------------
  control_exp_med <- sapply(unique(tab2p$cancer), FUN = function(cancer, tab2p) {
    tab2p_can <- tab2p[tab2p$cancer == cancer & tab2p$Classification_complex == "WT",]
    control_exp_med <- median(tab2p_can$exp_value, na.rm = T)
    return(control_exp_med)
  }, tab2p = tab2p)
  control_exp_meds <- data.frame(control_exp_med = control_exp_med, cancer = names(control_exp_med))
  tab2p <- merge(tab2p, control_exp_meds, by = c("cancer"), all.x = T)
  ## get the high quantile of expression values
  exp_high <- sapply(unique(tab2p$cancer), FUN = function(cancer, tab2p) {
    tab2p_can <- tab2p[tab2p$cancer == cancer,]
    exp_high <- quantile(x = tab2p_can$exp_value, probs = 0.9, na.rm = T)
    return(exp_high)
  }, tab2p = tab2p)
  names(exp_high) <-  unique(tab2p$cancer)
  exp_highs <- data.frame(exp_high = exp_high, cancer = unique(tab2p$cancer))
  tab2p <- merge(tab2p, exp_highs, by = c("cancer"), all.x = T)
  ## get the low quantile of expression values
  exp_low <- sapply(unique(tab2p$cancer), FUN = function(cancer, tab2p) {
    tab2p_can <- tab2p[tab2p$cancer == cancer,]
    exp_low <- quantile(x = tab2p_can$exp_value, probs = 0.1, na.rm = T)
    return(exp_low)
  }, tab2p = tab2p)
  names(exp_low) <-  unique(tab2p$cancer)
  exp_lows <- data.frame(exp_low = exp_low, cancer = unique(tab2p$cancer))
  tab2p <- merge(tab2p, exp_lows, by = c("cancer"), all.x = T)
  
  
  # add the position of missense --------------------------------------------
  tab2p$point_fill <- "grey50"
  tab2p$point_fill[tab2p$position_y > 0] <- "orange"
  tab2p$point_fill[(as.vector(tab2p$x) > as.vector(tab2p$exp_high)) & (tab2p$position_y > 0)] <- "red"
  tab2p$point_fill[grepl(x = tab2p$Classification_complex, pattern = "Truncation") & (tab2p$position_y < 0)] <- "blue"
  tab2p$point_fill[tab2p$CNA.TP53 == "deletion" & (tab2p$position_y < 0)] <- "lightblue"
  
  tab2p$cancer <- factor(tab2p$cancer, levels = c("LIHC", "BRCA", "OV", "CO", "UCEC", "CCRCC"))
  
  p = ggplot()
  p = p + geom_point(data = tab2p, mapping = aes(x=x, y=y, color = point_fill, shape = subtype), stroke = 1, alpha = 0.6, size = 2.5)
  p <- p + geom_hline(yintercept = c(hotspot - position_range[1]), linetype = 2, color = "grey50", alpha = 0.5)
  p = p + geom_vline(data = tab2p, mapping = aes(xintercept = control_exp_med), linetype = 2, color = "grey50", alpha = 0.5)
  p = p + geom_vline(data = tab2p, mapping = aes(xintercept = exp_high), linetype = 2, color = set1[1], alpha = 0.5)
  p = p + geom_vline(data = tab2p, mapping = aes(xintercept = exp_low), linetype = 2, color = set1[2], alpha = 0.5)
  p = p + geom_text_repel(data = tab2p[(as.vector(tab2p$x) > as.vector(tab2p$exp_high)) & tab2p$position %in% hotspot_proximal,], 
                          mapping = aes(x=x, y=y, label= as.character(HGVSp_Short)), color = set1[1], segment.color = "black",
                          force = 4, segment.size = 1, segment.alpha = 0.6, size = 4, alpha = 0.8)
  p = p + geom_text_repel(data = tab2p[(as.vector(tab2p$x) < as.vector(tab2p$exp_low)) & tab2p$position %in% hotspot_proximal,], 
                          mapping = aes(x=x, y=y, label= as.character(HGVSp_Short)), color = set1[2], segment.color = "black",
                          force = 4, segment.size = 1, segment.alpha = 0.6, size = 4, alpha = 0.8)
  p <- p + scale_y_continuous(breaks = c(c(hotspot - position_range[1]), -20, -40, -60), labels = c(as.character(hotspot), "Truncation", "Deletion", "WT"))
  p = p + scale_color_manual(values = c("red" = set1[1], "blue" = set1[2], "grey50" = "grey50", "orange" = "#FF7F00", "lightblue" = "#A6CEE3"))
  p = p + facet_grid(.~cancer, scales = "free", space = "free_y")
  p = p + theme_nogrid()
  p = p + labs(x = paste0(geneB, "_", phosphosite, " abundance(log2)"), y = paste0(geneA, " mutation position"))
  p = p + theme(axis.title.x = element_text(size = 15, face = "bold"),
                axis.title.y = element_text(size = 15, face = "bold"),
                axis.text.y = element_text(colour="black", size=8))
  p = p + theme(title = element_text(size = 18, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 8))
  p <- p + theme(panel.spacing.x = unit(0, "lines"),  panel.spacing.y = unit(0, "lines"))
  p <- p + theme(legend.position="bottom", legend.text = element_text(size = 5), legend.key.width = unit(1, "lines"))
  p <- p + guides(color = FALSE)
  p
  ggsave(file=fn, width = 5, height = 4, useDingbats = F)
  dev.off()
  
}

sup_tab %>%
  tail()



