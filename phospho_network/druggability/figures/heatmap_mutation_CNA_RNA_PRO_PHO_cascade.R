# Yige Wu @ WashU 2019 Jan
## plot a heatmap with genomics data and proteomics data and kinase-substrate score status for given pairs



# source ------------------------------------------------------------------
setwd(dir = "~/Box Sync/")
source('./cptac2p_analysis/preprocess_files/preprocess_files_shared.R')
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
library(dplyr)

plotpheatmap <- function(mat_value, color.palette, col_anno, ann_colors, width, height, mat_label, gaps_row_vector) {
  result <- tryCatch({
    if (is.null(ann_colors)) {
      my_heatmap <- pheatmap(mat_value, 
                             color = color.palette,
                             scale = "row",
                             na_col = "white", 
                             display_numbers = mat_label,
                             cluster_rows=F, cluster_cols=T, 
                             gaps_row = gaps_row_vector,
                             number_color = "black", fontsize_number = 15,
                             show_colnames = F, annotation_colors = ann_colors)
    } else {
      my_heatmap <- pheatmap(mat_value, 
                             color = color.palette,
                             annotation_col = col_anno,
                             scale = "row",
                             na_col = "white", 
                             display_numbers = mat_label,
                             cluster_rows=F, cluster_cols=T, 
                             gaps_row = gaps_row_vector,
                             number_color = "black", fontsize_number = 15,
                             show_colnames = F, annotation_colors = ann_colors)
    }

  }, error = function(err) {
    print(print(paste("MY_ERROR:  ",err)))
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row",
                           na_col = "white", 
                           display_numbers = mat_label,
                           cluster_rows=F, cluster_cols=F, show_colnames = F, annotation_colors = ann_colors)
    return(NA)
  }, warning = function(war) {
    print(print(paste("MY_WARNING:  ", war)))
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row",
                           na_col = "white", 
                           display_numbers = mat_label,
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
# color.palette <- rainbow(8)
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")

# input outlier results -----------------------------------------------
esscore_tab_outlier_drug <- fread(input = "./cptac2p/analysis_results/phospho_network/druggability/figures/grid_percent_patient_with_druggable_outlier/esscore_tab_outlier_drug_genes.txt", data.table = F)
pairs2plot <- list();
for (cancer in cancers2process) {
  pairs2plot[[cancer]] <- list()
  genes2process <- unique(esscore_tab_outlier_drug$GENE[esscore_tab_outlier_drug$cancer == cancer])
  genes2process <- intersect(oncogenes, genes2process)
  for (geneA in genes2process) {
    pairs2plot[[cancer]][[geneA]] <- unique(esscore_tab_outlier_drug$pair[esscore_tab_outlier_drug$GENE == geneA])
  }
}
# geneA <- "ERBB2"; pairs2plot[[geneA]] <- unique(esscore_tab_outlier_drug$pair[esscore_tab_outlier_drug$GENE == geneA])
# geneA <- "AKT1"; pairs2plot[[geneA]] <- unique(esscore_tab_outlier_drug$pair[esscore_tab_outlier_drug$GENE == geneA])
# geneA <- "MTOR"; pairs2plot[[geneA]] <- unique(esscore_tab_outlier_drug$pair[esscore_tab_outlier_drug$GENE == geneA])
# geneA <- "BRAF"; pairs2plot[[geneA]] <- unique(esscore_tab_outlier_drug$pair[esscore_tab_outlier_drug$GENE == geneA])
# geneA <- "MAP2K1"; pairs2plot[[geneA]] <- unique(esscore_tab_outlier_drug$pair[esscore_tab_outlier_drug$GENE == geneA])
# geneA <- "MAPK1"; pairs2plot[[geneA]] <- unique(esscore_tab_outlier_drug$pair[esscore_tab_outlier_drug$GENE == geneA])
genes_altered <- c("MAP3K1", "PIK3CA")
# bussiness ------------------------------------------------------------------
# for (cancer in c("BRCA")) {
for (cancer in c("BRCA", "UCEC", "CO", "CCRCC", "OV")) {
  subdir1 <- paste0(makeOutDir(resultD = resultD), paste(cancer, sep = "_"), "/")
  dir.create(subdir1)
  fn <- paste0(subdir1, cancer, "_outlier", outlier_sd, "SD_",  "withlegend.pdf")
  
  if (!file.exists(fn)) {
    mat_values_sup <- NULL
    mat_labels_sup <- NULL
    gaps_row_vector <- NULL
    ann_colors <- list()
    col_anno <- NULL
    
    
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
    maf <- loadMaf(cancer = cancer, maf_files = maf_files)
    mut_mat <- generate_somatic_mutation_matrix(pair_tab = names(pairs2plot[[cancer]]), maf = maf)
    cna_tab <- loadCNAstatus(cancer = cancer)
    partIDs <- colnames(pho_tab)[!(colnames(pho_tab) %in% c("Gene", "Phosphosite", "Peptide_ID"))]
    
    # file2input <- paste0("./cptac2p/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/kinase/esscore_tab_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt")
    file2input <- paste0("./cptac2p/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/kinase/esscore_tab_outlier", outlier_sd, "SD_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt")
    if (file.exists(file2input)) {
      ## input kinase-substrate score table
      ## input kinase-substrate score outlier status table
      esscore_tab_outlier <- fread(input = file2input, data.table = F)
      
      escore_tab <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/kinase/escore_tab_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
      escore_tab_scaled <- escore_tab[, partIDs]
      escore_tab_scaled <- scale_by_row(escore_tab_scaled)
      escore_tab_outlier <- cbind(esscore_tab_outlier[, c("pair", "SELF", "GENE")], (escore_tab_scaled > outlier_sd))
      # colnames(escore_tab_outlier) <- c("pair", "SELF", partIDs)
    } else {
      esscore_tab <- NULL
      esscore_tab_outlier <- NULL
      escore_tab <- NULL
      sscore_tab <- NULL
    }
    
    if (is.null(col_anno)) {
      col_anno <- data.frame(partID = partIDs)
    }
    
    if (cancer %in% c("BRCA", "CO", "UCEC")) {
      if (cancer == "BRCA") {
        subtypes <- partID2pam50(patientID_vector = partIDs, pam50_map = loadPAM50Map())
      } else if (cancer == "CO") {
        subtypes <- partID2MSI(patientID_vector = partIDs, subtype_map = loadMSIMap())
      } else if (cancer == "UCEC") {
        subtypes <- partID2UCECsubtype(patientID_vector = partIDs)
      }
      subtypes2merge <- data.frame(partID = partIDs, subtype = subtypes)
    }
    
    
    
    for (geneA in names(pairs2plot[[cancer]])) {
      pairs <- pairs2plot[[cancer]][[geneA]]
      pairs_mat <- str_split_fixed(string = pairs, pattern = ":", n = 3)
      site_ids <- paste0(pairs_mat[, 2], "_", pairs_mat[, 3])
      
      # make the matrix of values showing in heatmap ----------------------------
      df_value <- NULL
      
      pro_tab_tmp <- pro_tab[pro_tab$Gene %in% c(geneA),]
      if (nrow(pro_tab_tmp) > 0) {
        pro_tab.m <- melt(pro_tab_tmp, id.vars = "Gene")
        pro_tab.m %>% head()
        colnames(pro_tab.m) <- c("Gene", "partID", "exp_value")
        pro_tab.m$Phosphosite <- "PRO"
        pro_tab.m$row_id <- paste0(pro_tab.m$Gene, "_", "PRO")
        df_value <- rbind(df_value, pro_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value", "row_id")])
      }
      
      pho_tab_tmp <- pho_tab %>%
        mutate(site_id = paste0(Gene, "_", Phosphosite)) %>%
        filter(site_id %in% site_ids)
      if (nrow(pho_tab_tmp) > 0) {
        pho_tab_tmp <- pho_tab_tmp[, c("Gene", "Phosphosite", partIDs)]
        pho_tab.m <- melt(pho_tab_tmp, id.vars = c("Gene", "Phosphosite"))
        pho_tab.m %>% head()
        colnames(pho_tab.m) <- c("Gene", "Phosphosite", "partID", "exp_value")
        pho_tab.m$row_id <- paste0(pho_tab.m$Gene, "_", pho_tab.m$Phosphosite)
        df_value <- rbind(df_value, pho_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value", "row_id")])
      }
      
      phog_tab_tmp <- phog_tab[phog_tab$Gene %in% c(geneA),]
      if (!is.null(phog_tab_tmp)) {
        if (nrow(phog_tab) > 0) {
          phog_tab.m <- melt(phog_tab_tmp, id.vars = "Gene")
          phog_tab.m %>% head()
          colnames(phog_tab.m) <- c("Gene", "partID", "exp_value")
          phog_tab.m$Phosphosite <- "collapsed_PHO"
          phog_tab.m$row_id <- paste0(phog_tab.m$Gene, "_", phog_tab.m$Phosphosite)
          df_value <- rbind(df_value, phog_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value", "row_id")])
        }
      }
      
      ## if it's outlier in k-s score
      esscore_tab_outlier_tmp <- esscore_tab_outlier[esscore_tab_outlier$pair %in% pairs,]
      
      if (!is.null(esscore_tab_outlier_tmp)) {
        if (nrow(esscore_tab_outlier_tmp) > 0) {
          esscore_tab_outlier.m <- melt(esscore_tab_outlier_tmp[, c("pair", "SUB_GENE", "SUB_MOD_RSD", partIDs)], 
                                        id.vars = c("pair", "SUB_GENE", "SUB_MOD_RSD"))
          esscore_tab_outlier.m <- esscore_tab_outlier.m %>%
            filter(!is.na(value) & value == T) %>%
            mutate(site_id = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
            select(site_id, variable, value)
          if (nrow(esscore_tab_outlier.m) > 0) {
            colnames(esscore_tab_outlier.m) <- c("site_id", "partID", "esscore_outlier")
            df_value <- merge(df_value, esscore_tab_outlier.m[, c("partID", "site_id", "esscore_outlier")], 
                              by.x = c("partID", "row_id"), 
                              by.y = c("partID", "site_id"),
                              all.x = T)
            cis_outlier_partIDs <- df_value$partID[df_value$Gene == geneA & df_value$Phosphosite != "PRO" & !is.na(df_value$esscore_outlier)  & (df_value$esscore_outlier == T)]
            df_value$esscore_outlier[df_value$Gene == geneA & df_value$Phosphosite == "PRO" & (df_value$partID %in% cis_outlier_partIDs)] <- T
            
            trans_outlier_partIDs <- df_value$partID[df_value$Gene != geneA & df_value$Phosphosite != "PRO" & !is.na(df_value$esscore_outlier)  & (df_value$esscore_outlier == T)]
            df_value$esscore_outlier[df_value$Gene == geneA & df_value$Phosphosite == "collapsed_PHO" & (df_value$partID %in% trans_outlier_partIDs)] <- T
            # col_anno$is_ks_outlier <- as.character(rownames(col_anno) %in% df_value$partID[df_value$esscore_outlier])
          }
        }
      }
      
      df_value$exp_value <- as.numeric(as.vector(df_value$exp_value))
      df_value <- unique(df_value)
      df_value <- df_value %>%
        filter(!is.na(df_value$exp_value)) %>%
        mutate(id = paste0(row_id, ":", partID)) %>%
        filter(!duplicated(id)) %>%
        unique
      
      ## make the matrix for the heatmap body
      mat_value <- dcast(data = df_value, row_id ~ partID, value.var = "exp_value")
      
      mat_value %>% head()
      rownames(mat_value) <- mat_value$row_id
      mat_value <- as.matrix(mat_value[,-1])
      head(mat_value)
      # mat_value <- mat_value[rowSums(!is.na(mat_value)) > 0.5*ncol(mat_value),]
      
      mat_label <- dcast(data = df_value, row_id ~ partID, value.var = "esscore_outlier")
      rownames(mat_label) <- mat_label$row_id
      mat_label <- as.matrix(mat_label[,-1])
      
      top_rows <- sort(rowSums(!is.na(mat_label)), decreasing = T)
      top_rows <- names(top_rows)[1:min(length(top_rows), 10)]
      top_rows <- c(top_rows[!grepl(x = top_rows, pattern = geneA)], top_rows[grepl(x = top_rows, pattern = geneA)])
      
      mat_value <- mat_value[intersect(top_rows, rownames(mat_value)),]
      
      # mat_value <- mat_value[,colSums(!is.na(mat_value)) > 0]
      mat_value %>% head()
      
      # make the matrix of column annotation showing in heatmap ----------------------------
      
      ## if it's outlier in k-s score
      esscore_tab_outlier_tmp <- esscore_tab_outlier[esscore_tab_outlier$pair %in% pairs,]
      if (!is.null(esscore_tab_outlier_tmp)) {
        if (nrow(esscore_tab_outlier_tmp) > 0) {
          esscore_tab_outlier.m <- melt(esscore_tab_outlier_tmp[, c("pair", "SUB_GENE", "SUB_MOD_RSD", partIDs)], 
                                        id.vars = c("pair", "SUB_GENE", "SUB_MOD_RSD"))
          esscore_tab_outlier.m <- esscore_tab_outlier.m %>%
            filter(!is.na(value) & value == T) %>%
            mutate(site_id = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
            filter(site_id %in% rownames(mat_value)) %>%
            select(variable, value, pair, site_id)
          
          site_ids_ordered <- sort(unique(esscore_tab_outlier.m$site_id))
          site_ids_ordered <- c(site_ids_ordered[!grepl(x = site_ids_ordered, pattern = geneA)], site_ids_ordered[grepl(x = site_ids_ordered, pattern = geneA)])
          esscore_tab_outlier.m$site_id <- factor(esscore_tab_outlier.m$site_id, levels = site_ids_ordered)
          esscore_tab_outlier.m <- esscore_tab_outlier.m[order(esscore_tab_outlier.m$site_id),]
          
          if (nrow(esscore_tab_outlier.m) > 0) {
            for (pair_tmp in unique(esscore_tab_outlier.m$pair)) {
              esscore_tab_outlier2merge <- esscore_tab_outlier.m %>%
                filter(pair == pair_tmp) %>%
                select(variable, value)
              colnames(esscore_tab_outlier2merge) <- c("partID", pair_tmp)
              col_anno <- merge(col_anno, esscore_tab_outlier2merge, by = c("partID"), all.x = T)
              col_anno[, pair_tmp] <- ifelse(is.na(col_anno[, pair_tmp]), "FALSE", "TRUE")
              ann_colors[[pair_tmp]] <- c("TRUE" = "#ce1256", "FALSE" = "#d9d9d9")
            }
          }
        }
      }
      
      phog_outlier_tmp <- escore_tab_outlier[escore_tab_outlier$SELF == "trans" & escore_tab_outlier$GENE == geneA,]
      phog_outlier_tmp <- phog_outlier_tmp[1,]
      if (!is.null(phog_outlier_tmp)) {
        if (nrow(phog_outlier_tmp) > 0) {
          phog_outlier.m <- melt(phog_outlier_tmp, 
                                 id.vars = colnames(phog_outlier_tmp)[!(colnames(phog_outlier_tmp) %in% partIDs)])
          phog_outlier2merge <- phog_outlier.m %>%
            select(variable, value)
          colnames(phog_outlier2merge) <- c("partID", paste0("PHO.", geneA))
          col_anno <- merge(col_anno, phog_outlier2merge, by = c("partID"), all.x = T)
          col_anno[, paste0("PHO.", geneA)] <- ifelse(col_anno[, paste0("PHO.", geneA)] == TRUE & !is.na(col_anno[, paste0("PHO.", geneA)]), "TRUE", "FALSE")
          ann_colors[[paste0("PHO.", geneA)]] <- c("TRUE" = "#e7298a", "FALSE" = "#d9d9d9")
        }
      }
      
      pro_outlier_tmp <- escore_tab_outlier[escore_tab_outlier$SELF == "cis" & escore_tab_outlier$GENE == geneA,]
      pro_outlier_tmp <- pro_outlier_tmp[1,]
      if (!is.null(pro_outlier_tmp)) {
        if (nrow(pro_outlier_tmp) > 0) {
          pro_outlier.m <- melt(pro_outlier_tmp, 
                                id.vars = colnames(pro_outlier_tmp)[!(colnames(pro_outlier_tmp) %in% partIDs)])
          pro_outlier2merge <- pro_outlier.m %>%
            select(variable, value)
          colnames(pro_outlier2merge) <- c("partID", paste0("PRO.", geneA))
          col_anno <- merge(col_anno, pro_outlier2merge, by = c("partID"), all.x = T)
          col_anno[, paste0("PRO.", geneA)] <- ifelse(col_anno[, paste0("PRO.", geneA)] == TRUE & !is.na(col_anno[, paste0("PRO.", geneA)]), "TRUE", "FALSE")
          ann_colors[[paste0("PRO.", geneA)]] <- c("TRUE" = "#df65b0", "FALSE" = "#d9d9d9")
        }
      }
      
      ## CNA needs to show both geneA and geneB; levels: amplification, deletion, neutral
      cna_tab_tmp <- cna_tab[cna_tab$gene %in% geneA, ]
      if (nrow(cna_tab_tmp) > 0) {
        cna_tab.m <- melt(cna_tab_tmp, id.vars = "gene")
        colnames(cna_tab.m) <- c("gene", "partID", "CNA")
        
        for (gene in unique(cna_tab.m$gene[cna_tab.m$CNA != "neutral"])) {
          cna_mat2merge <- cna_tab.m[cna_tab.m$gene == gene, c("partID", "CNA")]
          colnames(cna_mat2merge) <- c("partID", paste0("CNA.", gene))
          col_anno <- merge(col_anno, cna_mat2merge, by = c("partID"), all.x = T)
        }
      } else {
        print("no CNA!")
      }
      
      mut_mat_tmp <- mut_mat[mut_mat$Hugo_Symbol ==  geneA,]
      if (nrow(mut_mat_tmp) > 0){
        mut_mat.m <- melt(mut_mat_tmp, id.vars = "Hugo_Symbol")
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
          col_anno[is.na(col_anno[, paste0("mutation.", gene)]), paste0("mutation.", gene)] <- "wild_type"
        }
      } else {
        print("no mutation!")
        col_anno[,paste0("mutation.", gene)] <- "wild_type"
      }
      
      if (paste0("CNA.", geneA) %in% colnames(col_anno)) {
        ann_colors[[paste0("CNA.", geneA)]] <-  c(amplification = "#c994c7", deletion = "#d9d9d9", "neutral" = "#d9d9d9")
      }
      
      if (paste0("mutation.", geneA) %in% colnames(col_anno)) {
        ann_colors[[paste0("mutation.", geneA)]] <- c(missense = "#fa9fb5", truncation = "#377EB8", wild_type = "#d9d9d9", "missense&truncation" = "#6A3D9A", other_mutation = "#FF7F00", silent = "#33A02C")
      }
      
      col_anno %>% head()
      rownames(col_anno) <- col_anno$partID
      
      
      ## order the matrix column
      partIDs_ordered <- intersect(as.vector(rownames(col_anno)), colnames(mat_value))
      # partIDs_ordered <- names(sort(colSums(mat_label[, partIDs_ordered], na.rm = T)))
      col_anno <- col_anno[partIDs_ordered,]
      mat_value <- mat_value[, partIDs_ordered]
      
      ## order the matrix rows
      # row_order <- c(paste0(geneA, "_", c("RNA", "PRO", "collapsed_PHO")), site_ids)
      row_order <- c(paste0(geneA, "_", "PRO"), site_ids[grepl(pattern = geneA, x = site_ids)], 
                     paste0(geneA, "_", "collapsed_PHO"), site_ids[!grepl(pattern = geneA, x = site_ids)])
      mat_value <- mat_value[intersect(row_order, rownames(mat_value)),]
      
      mat_label <- mat_label[rownames(mat_value), colnames(mat_value)]
      mat_label[!is.na(mat_label)] <- "x"
      mat_label[is.na(mat_label)] <- ""
      
      # plotting ----------------------------------------------------------------
      
      ## reformat the outlier status column
      if ("esscore_outlier" %in% colnames(col_anno)) {
        tmp <- vector(mode = "character", length = nrow(col_anno))
        tmp[is.na(col_anno$esscore_outlier)] <- NA
        tmp[!is.na(col_anno$esscore_outlier) & col_anno$esscore_outlier == TRUE] <- "TRUE"
        tmp[!is.na(col_anno$esscore_outlier) & col_anno$esscore_outlier == FALSE] <- "FALSE"
        col_anno$esscore_outlier <- tmp
      }
      mat_values_sup <- rbind(mat_values_sup, mat_value)
      mat_labels_sup <- rbind(mat_labels_sup, mat_label)
      if (is.null(gaps_row_vector)) {
        gaps_row_vector <- c(gaps_row_vector, nrow(mat_value))
      } else {
        if (any(!grepl(rownames(mat_value), pattern = geneA))) {
          gaps_row_vector <- c(gaps_row_vector, 
                               (gaps_row_vector[length(gaps_row_vector)] + min(which(!grepl(rownames(mat_value), pattern = geneA))) - 2),
                               (nrow(mat_value) + gaps_row_vector[length(gaps_row_vector)]))
          
        } else {
          gaps_row_vector <- c(gaps_row_vector, 
                               (nrow(mat_value) + gaps_row_vector[length(gaps_row_vector)]))
        }
      }
    }
    
    ## subtype info, need to adapt when subtype info is absent
    if (!("subtype" %in% colnames(col_anno))) {
      col_anno <- merge(col_anno, subtypes2merge, by = c("partID"), all.x = T)
    }
    
    if ("subtype" %in% colnames(col_anno)) {
      subtype_colors <- set2[1:length(unique(subtypes))]
      names(subtype_colors) <- unique(subtypes)
      ann_colors[["subtype"]] <- subtype_colors
    }
    
    rownames(col_anno) <- col_anno$partID
    col_anno$partID <- NULL
    
    if (!("subtype" %in% colnames(col_anno))) {
      plotpheatmap(mat_value = mat_values_sup, color.palette = color.palette, ann_colors = ann_colors, width = 20, height = 0.1*nrow(mat_values_sup) + 3, mat_label = mat_labels_sup, gaps_row_vector = gaps_row_vector, col_anno = NULL)
      
    } else {
      plotpheatmap(mat_value = mat_values_sup, color.palette = color.palette, col_anno = col_anno, ann_colors = ann_colors, 20, 20, mat_labels_sup, gaps_row_vector)
    }
  }
}


