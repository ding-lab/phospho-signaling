# Yige Wu @ WashU 2019 Jan
## plot a heatmap with genomics data and proteomics data and kinase-substrate score status for given pairs



# source ------------------------------------------------------------------
setwd(dir = "~/Box Sync/")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
library(dplyr)

plotpheatmap <- function(mat_value, color.palette, col_anno, ann_colors, width, height) {
  result <- tryCatch({
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row",
                           na_col = "white", 
                           cluster_rows=F, cluster_cols=T, show_colnames = F, annotation_colors = ann_colors)
  }, error = function(err) {
    print(print(paste("MY_ERROR:  ",err)))
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row",
                           na_col = "white", 
                           cluster_rows=F, cluster_cols=F, show_colnames = F, annotation_colors = ann_colors)
    return(NA)
  }, warning = function(war) {
    print(print(paste("MY_WARNING:  ", war)))
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row",
                           na_col = "white", 
                           cluster_rows=F, cluster_cols=F, show_colnames = F, annotation_colors = ann_colors)
    
    return(NULL)
  }, finally = {
    save_pheatmap_pdf(x = my_heatmap, 
                      filename = fn, 
                      width = width, height = height)
  })
  return(result)
}

# set variables -----------------------------------------------------------
## variables for inputting outlier score table
enzyme_type <- "kinase"
reg_nonNA <- 20
outlier_sd <- 1.5

## plotting paramters
cap <- 3
breaks = seq(-(cap),cap, by=0.2)
## add color palette
color.palette <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(length(breaks))


# input fisher test results -----------------------------------------------
fisher_stat_cancers <- fread(input = "./cptac2p/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/fisher_stat_tab_outlier1.5SD_all_cancers_kinase_reg_nonNA20.txt", data.table = F)

# Filter by gene list -----------------------------------------------------
# genes2plot <- unique(c(unlist(SMGs), "ESR1", "MYC", "EGFR", "ERBB2"))
# fisher_stat_cancers2plot <- fisher_stat_cancers[fisher_stat_cancers$p.value < 0.05 & (fisher_stat_cancers$GENE %in% genes2plot | fisher_stat_cancers$SUB_GENE %in% genes2plot),]
fisher_stat_cancers2plot <- fisher_stat_cancers[fisher_stat_cancers$fdr < 0.1 & fisher_stat_cancers$subtype == "LumB",]
fisher_stat_cancers2plot <- fisher_stat_cancers[fisher_stat_cancers$fdr < 0.1 & fisher_stat_cancers$subtype == "LumA",]

fisher_stat_cancers2plot <- fisher_stat_cancers %>%
  filter(subtype == "POLE", p.value < 0.05, num_is.outlier_is.subtype >= 4)

fisher_stat_cancers2plot <- fisher_stat_cancers %>%
  filter(subtype == "Serous", fdr < 0.05, num_is.outlier_is.subtype >= 4)

fisher_stat_cancers2plot <- fisher_stat_cancers %>%
  filter(subtype == "Endometrioid", p.value < 0.05, num_is.outlier_is.subtype >= 4)

fisher_stat_cancers2plot <- fisher_stat_cancers[fisher_stat_cancers$SUB_GENE == "ERF",]


pairs2plot <- unique(fisher_stat_cancers2plot$pair)
# pairs2plot <- c("TP53:ESR1:PRO:MDM2")
# pairs2plot <- c("AKT1:BAD:S99:PIK3CA:PTEN:PIK3R1", "AKT1:GSK3B:S9:PIK3CA:PTEN", "PTEN:GSK3B:S9:PIK3CA:AKT1")
# pairs2plot <- c("ERBB2:ERBB2:Y1109", "ERBB2:ERBB2:S1083", "CDK1:TP53:S315", "DYRK1A:CCND1:T286", "MAPK9:JUN:S63")

mmr_gene <- read_delim("~/Box Sync/Ding_Lab/Projects_Current/MSI/MSI_shared_data/mmr_gene.txt","\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
mmr_gene <- as.vector(mmr_gene$X1)
pairs2plot <- c(paste0("TP53:MSH2:", paste0(mmr_gene, collapse = ":")))

pairs2plot <- c("PRKD1:CTNNB1:S552")

# pairs2plot <- c("MAPK1:MAPT:S721:MAP3K1:PIK3CA")
# pairs2plot <- c("MAPK1:MAPT:S721:ESR1")
pairs2plot <- c("GATA3:EGFR:PRO")
pairs2plot <- c("PIK3CA:AKT2:Y38")
pairs2plot <- c("TP53:TP53BP1:S1623:PRO:PLK1:CHEK1:CHEK2:MET:CDK1")
pairs2plot <- c("TP53:CDK1:PRO")
pairs2plot <- c("TP53:RB1:S807:CDK1:CDK4")
pairs2plot <- c("TP53:PLK1:PRO")
pairs2plot <- c("TP53:CHEK1:PRO")
pairs2plot <- c("TP53:CHEK2:PRO")
pairs2plot <- c("TP53:TP53BP1:S1623:PLK1")
pairs2plot <- c("TP53:TP53BP1:S1623:CDK1")
pairs2plot <- c("TP53:CTNNB1:PRO")

# bussiness ------------------------------------------------------------------
for (pair in pairs2plot) {
  geneA <- str_split(string = pair, pattern = ":")[[1]][1]
  geneB <- str_split(string = pair, pattern = ":")[[1]][2]
  phosphosite <- str_split(string = pair, pattern = ":")[[1]][3]
  ## geneC would be annotational mutations to show
  num_genes <- length(str_split(string = pair, pattern = ":")[[1]])
  geneC <- str_split(string = pair, pattern = ":")[[1]][4:num_genes]
  
  subdir1 <- paste0(makeOutDir(resultD = resultD), paste(geneA, geneB, phosphosite, sep = "_"), "/")
  dir.create(subdir1)
  for (cancer in c("UCEC", "BRCA", "CCRCC", "CO", "OV")) {
    fn <- paste0(subdir1, paste(geneA, geneB, phosphosite, sep = "_"), "_", cancer, "_outlier", outlier_sd, "SD_",  "withlegend.pdf")
    
    if (!file.exists(fn)) {
      ann_colors <- list()
      
      # input data first because different for each cancer type --------------------------------------------------------------
      ## input mutation matrix
      mut_mat <- fread(paste0("./cptac2p/analysis_results/phospho_network/genoalt/tables/test_mut_impact_proteome/", cancer, "_somatic_mutation.txt"), data.table = F)
      ## mutation needs to show both geneA and geneB
      # mut_mat <- mut_mat[mut_mat$Hugo_Symbol %in% c(geneA, geneB, geneC),]
      mut_mat <- mut_mat[mut_mat$Hugo_Symbol %in% c(geneA, geneB),]
      
      ## input CNA matrix
      cna_tab <- loadCNAstatus(cancer = cancer)
      # cna_tab <- cna_tab[cna_tab$gene %in% c(geneA, geneB, geneC), ]
      cna_tab <- cna_tab[cna_tab$gene %in% c(geneA, geneB), ]
      
      file2input <- paste0("./cptac2p/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/kinase/esscore_tab_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt")
      if (file.exists(file2input)) {
        ## input kinase-substrate score table
        esscore_tab <- fread(input = file2input, data.table = F)
        esscore_tab <- esscore_tab[esscore_tab$pair == pair,]
        
        ## input kinase-substrate score outlier status table
        esscore_tab_outlier <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/kinase/esscore_tab_outlier", outlier_sd, "SD_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
        esscore_tab_outlier <- esscore_tab_outlier[esscore_tab_outlier$pair == pair,]
        
        ## input kinase score table
        escore_tab <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/kinase/escore_tab_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
        escore_tab <- escore_tab[escore_tab$pair == pair,]
        
        ## input substrate score table
        sscore_tab <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/kinase/sscore_tab_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
        sscore_tab <- sscore_tab[sscore_tab$pair == pair,]
      } else {
        esscore_tab <- NULL
        esscore_tab_outlier <- NULL
        escore_tab <- NULL
        sscore_tab <- NULL
      }
      
      ## load RNA
      rna_tab <- loadRNA(cancer = cancer)
      rna_tab <- rna_tab[rna_tab$gene %in% c(geneA, geneB),]
      
      ## input protein data
      if (cancer %in% c("BRCA", "OV", "CO")) {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
        phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
        
      } else if (cancer == "UCEC") {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
        phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
        
      } else if (cancer == "CCRCC") {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
        phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
      } else if (cancer == "LIHC") {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
        phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
      }
      pro_tab <- pro_tab[pro_tab$Gene %in% c(geneA, geneB, geneC),]
      pho_tab <- pho_tab[pho_tab$Gene == geneB & pho_tab$Phosphosite == phosphosite,]
      phog_tab <- phog_tab[phog_tab$Gene %in% c(geneA, geneB),]
      
      # make the annotation columns for each sample -----------------------------
      partIDs <- colnames(pho_tab)[!(colnames(pho_tab) %in% c("Gene", "Phosphosite", "Peptide_ID"))]
      col_anno <- data.frame(partID = partIDs)
      
      if (nrow(mut_mat) > 0){
        mut_mat.m <- melt(mut_mat, id.vars = "Hugo_Symbol")
        mut_mat.m %>% head()
        colnames(mut_mat.m) <- c("Gene", "partID", "variant_class")
        
        ## distinguish by missense and truncation
        mut_mat.m$variant_class[is.na(mut_mat.m$variant_class)] <- ""
        mut_mat.m$variant_class_sim <- "other_mutation"
        mut_mat.m$variant_class_sim[mut_mat.m$variant_class == ""] <- "wild_type"
        mut_mat.m$variant_class_sim[mut_mat.m$variant_class  == "Silent"] <- "silent"
        mut_mat.m$variant_class_sim[grepl(x = mut_mat.m$variant_class, pattern = "Missense_Mutation")] <- "missense"
        mut_mat.m$variant_class_sim[grepl(x = mut_mat.m$variant_class, pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del")] <- "truncation"
        mut_mat.m$variant_class_sim[sapply(X = mut_mat.m$variant_class, FUN = function(v) (grepl(pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del", x = v) & grepl(pattern = "Missense_Mutation", x = v)))] <- "missense&truncation"
        
        for (gene in unique(mut_mat.m$Gene[mut_mat.m$variant_class_sim != "wild_type"])) {
          mut_mat2merge <- mut_mat.m[mut_mat.m$Gene == gene, c("partID", "variant_class_sim")]
          colnames(mut_mat2merge) <- c("partID", paste0("mutation.", gene))
          col_anno <- merge(col_anno, mut_mat2merge, by = c("partID"), all.x = T)
        }
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
      
      ## subtype info, need to adapt when subtype info is absent
      if (cancer %in% c("BRCA", "CO", "UCEC")) {
        if (cancer == "BRCA") {
          subtypes <- partID2pam50(patientID_vector = partIDs, pam50_map = loadPAM50Map())
        } else if (cancer == "CO") {
          subtypes <- partID2MSI(patientID_vector = partIDs, subtype_map = loadMSIMap())
        } else if (cancer == "UCEC") {
          subtypes <- partID2UCECsubtype(patientID_vector = partIDs)
        }
        subtypes2merge <- data.frame(partID = partIDs, subtype = subtypes)
        col_anno <- merge(col_anno, subtypes2merge, by = c("partID"), all.x = T)
        
      }
      
      ## if it's outlier in k-s score
      if (!is.null(esscore_tab)) {
        if (nrow(esscore_tab) > 0) {
          esscore_tab.m <- melt(esscore_tab[, c("pair", partIDs)], id.vars = "pair")
          esscore_tab.m <- esscore_tab.m[!is.na(esscore_tab.m$value),]
          if (nrow(esscore_tab.m) > 0) {
            colnames(esscore_tab.m) <- c("pair", "partID", "esscore")
            col_anno <- merge(col_anno, esscore_tab.m[, c("partID", "esscore")], by = c("partID"), all.x = T)
            rm(esscore_tab.m)
            
            esscore_tab_outlier.m <- melt(esscore_tab_outlier[, c("pair", partIDs)], id.vars = "pair")
            esscore_tab_outlier.m <- esscore_tab_outlier.m[!is.na(esscore_tab_outlier.m$value),]
            colnames(esscore_tab_outlier.m) <- c("pair", "partID", "esscore_outlier")
            col_anno <- merge(col_anno, esscore_tab_outlier.m[, c("partID", "esscore_outlier")], by = c("partID"), all.x = T)
            rm(esscore_tab_outlier.m)
            
            ## add kinase score outlier status
            escore_tab.m <- melt(escore_tab[, c("pair", partIDs)], id.vars = "pair")
            escore_tab.m <- escore_tab.m[!is.na(escore_tab.m$value),]
            colnames(escore_tab.m) <- c("pair", "partID", "escore")
            escore_tab.m$escore_outlier <- ifelse(escore_tab.m$escore > outlier_sd, "TRUE", "FALSE")
            col_anno <- merge(col_anno, escore_tab.m[, c("partID", "escore_outlier")], by = c("partID"), all.x = T)
            
            ## add substrate score outlier status
            sscore_tab.m <- melt(sscore_tab[, c("pair", partIDs)], id.vars = "pair")
            sscore_tab.m <- sscore_tab.m[!is.na(sscore_tab.m$value),]
            colnames(sscore_tab.m) <- c("pair", "partID", "sscore")
            sscore_tab.m$sscore_outlier <- ifelse(sscore_tab.m$sscore > outlier_sd, "TRUE", "FALSE")
            col_anno <- merge(col_anno, sscore_tab.m[, c("partID", "sscore_outlier")], by = c("partID"), all.x = T)
          } 
        }
        
      }
      
      ## order samples
      col_anno %>% head()
      
      for (gene in unique(c(geneA, geneB, geneC))) {
        if (paste0("CNA.", gene) %in% colnames(col_anno)) {
          col_anno <- col_anno[order(col_anno[, paste0("CNA.", gene)], decreasing = T),]
          ann_colors[[paste0("CNA.", gene)]] <-  c(amplification = "#E41A1C", deletion = "#377EB8", "neutral" = "grey")
        }
        if (paste0("mutation.", gene) %in% colnames(col_anno)) {
          col_anno <- col_anno[order(col_anno[, paste0("mutation.", gene)], decreasing = T),]
          ann_colors[[paste0("mutation.", gene)]] <- c(missense = "#E41A1C", truncation = "#377EB8", wild_type = "white", "missense&truncation" = "#6A3D9A", other_mutation = "#FF7F00", silent = "#33A02C")
        }
      }
      if ("subtype" %in% colnames(col_anno)) {
        col_anno <- col_anno[order(col_anno$subtype, decreasing = T),]
        subtype_colors <- set2[1:length(unique(subtypes))]
        names(subtype_colors) <- unique(subtypes)
        ann_colors[["subtype"]] <- subtype_colors
      }
      if ("esscore_outlier" %in% colnames(col_anno)) {
        ann_colors[["esscore_outlier"]] <- c("TRUE" = "#E41A1C", "FALSE" = "grey")
        ann_colors[["escore_outlier"]] <- c("TRUE" = "#E41A1C", "FALSE" = "grey")
        ann_colors[["sscore_outlier"]] <- c("TRUE" = "#E41A1C", "FALSE" = "grey")
        tmp <- color.palette; names(tmp) <- seq(from = -2, to = 2, length.out = length(color.palette))
        ann_colors[["esscore"]] <- tmp
        # col_anno <- col_anno[order(col_anno$esscore, decreasing = T),]
      }
      
      
      # col_anno <- col_anno[order(col_anno$variant_class_sim),]
      col_anno %>% head()
      rownames(col_anno) <- col_anno$partID
      
      # make the matrix of values showing in heatmap ----------------------------
      sup_tab_can <- NULL
      
      if (nrow(rna_tab) > 0) {
        rna_tab.m <- melt(rna_tab, id.vars = "gene")
        colnames(rna_tab.m) <- c("Gene", "partID", "exp_value")
        rna_tab.m$Phosphosite <- "RNA"
        sup_tab_can <- rbind(sup_tab_can, rna_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value")])
      }
      
      if (nrow(pro_tab) > 0) {
        pro_tab.m <- melt(pro_tab, id.vars = "Gene")
        pro_tab.m %>% head()
        colnames(pro_tab.m) <- c("Gene", "partID", "exp_value")
        pro_tab.m$Phosphosite <- "PRO"
        sup_tab_can <- rbind(sup_tab_can, pro_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value")])
      }
      
      if (nrow(pho_tab) > 0) {
        pho_tab <- pho_tab[!(colnames(pho_tab) %in% c("Peptide_ID"))]
        pho_tab.m <- melt(pho_tab, id.vars = c("Gene", "Phosphosite"))
        pho_tab.m %>% head()
        colnames(pho_tab.m) <- c("Gene", "Phosphosite", "partID", "exp_value")
        sup_tab_can <- rbind(sup_tab_can, pho_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value")])
      }
      
      # if (!is.null(phog_tab)) {
      #   if (nrow(phog_tab) > 0) {
      #     phog_tab.m <- melt(phog_tab, id.vars = "Gene")
      #     phog_tab.m %>% head()
      #     colnames(phog_tab.m) <- c("Gene", "partID", "exp_value")
      #     phog_tab.m$Phosphosite <- "collapsed_PHO"
      #     sup_tab_can <- rbind(sup_tab_can, phog_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value")])
      #   }
      #   
      # }
      
      sup_tab_can$id_row <- paste0(sup_tab_can$Gene, "_", sup_tab_can$Phosphosite)
      sup_tab_can$exp_value <- as.numeric(as.vector(sup_tab_can$exp_value))
      sup_tab_can <- unique(sup_tab_can)
      
      ## make the matrix for the heatmap body
      df_value <- dcast(data = sup_tab_can, id_row ~ partID, value.var = "exp_value")
      
      df_value %>% head()
      mat_value <- as.matrix(df_value[,-1])
      rownames(mat_value) <- df_value$id_row
      head(mat_value)
      # mat_value <- mat_value[,colSums(!is.na(mat_value)) > 0]
      mat_value %>% head()
      
      ## order the matrix column
      partIDs_ordered <- intersect(as.vector(rownames(col_anno)), colnames(mat_value))
      col_anno$partID <- NULL
      mat_value <- mat_value[, partIDs_ordered]
      
      ## order the matrix rows
      row_order <- c(paste0(rep(c(geneA, geneB, geneC), 3), "_", rep(c("RNA", "PRO", "collapsed_PHO"), c(2,2,2))), paste0(geneB, "_", phosphosite))
      mat_value <- mat_value[intersect(row_order, rownames(mat_value)),]
      
      # plotting ----------------------------------------------------------------
      
      ## reformat the outlier status column
      if ("esscore_outlier" %in% colnames(col_anno)) {
        tmp <- vector(mode = "character", length = nrow(col_anno))
        tmp[is.na(col_anno$esscore_outlier)] <- NA
        tmp[!is.na(col_anno$esscore_outlier) & col_anno$esscore_outlier == TRUE] <- "TRUE"
        tmp[!is.na(col_anno$esscore_outlier) & col_anno$esscore_outlier == FALSE] <- "FALSE"
        col_anno$esscore_outlier <- tmp
      }
      

      plotpheatmap(mat_value, color.palette, col_anno, ann_colors, 20, 10)

    }
  }
}


