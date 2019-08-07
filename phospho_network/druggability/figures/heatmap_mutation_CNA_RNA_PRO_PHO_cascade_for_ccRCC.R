# Yige Wu @ WashU 2019 Jan
## plot a heatmap with genomics data and proteomics data and kinase-substrate score status for given pairs

# source ------------------------------------------------------------------
setwd(dir = "~/Box Sync/")
source('./cptac2p_analysis/preprocess_files/preprocess_files_shared.R')
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
library(dplyr)
library(UpSetR)

# set variables -----------------------------------------------------------
## the list of druggable genes to output the sample annotation table
drug_genes2w <- c("MET", "BRAF", "MAP2K1", "MAPK1")
mut_genes2w <- SMGs[["CCRCC"]]

## variables for inputting outlier score table
enzyme_type <- "kinase"
reg_nonNA <- 20
outlier_sd <- 1.5
expressionType2colorValue <- function(vector_expressionType) {
  vector_colorValue <- vector(mode = "numeric", length = length(vector_expressionType))
  vector_colorValue[vector_expressionType == "KS"] <- 1
  vector_colorValue[vector_expressionType == "PHO"] <- 2
  vector_colorValue[vector_expressionType == "PRO"] <- 3
  vector_colorValue[vector_expressionType == "RNA"] <- 4
  vector_colorValue[vector_expressionType == "CNA"] <- 5
  vector_colorValue[vector_expressionType == "Mutation"] <- 6
  
  return(vector_colorValue)
}

## plotting paramters
cap <- 3
breaks = seq(-(cap),cap, by=0.2)
## add color palette
color.palette <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(length(breaks))
# color.palette <- rainbow(8)
# cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")
cancers2process <- c("CCRCC")


# input regression --------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)

regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression <- annotate_ks_source(regression = regression)

# input outlier results -----------------------------------------------
esscore_tab_outlier_drug <- fread(input = "./cptac2p/analysis_results/phospho_network/druggability/figures/grid_percent_patient_with_druggable_outlier/esscore_tab_outlier_drug_genes.txt", data.table = F)
esscore_tab_outlier_drug <- annotate_ks_source(regression = esscore_tab_outlier_drug)
esscore_tab_outlier_drug <- esscore_tab_outlier_drug %>%
  mutate(pair_cancer = paste0(pair, ":", cancer)) %>%
  mutate(regulated = (pair_cancer %in% regression$pair_cancer[regression$regulated == T]))

# input the druggable pairs to be annotated in the sample table-------------------------------------------------------------
cancer_tmp <- "CCRCC"
pairs2plot <- list()
pairs2plot[[cancer_tmp]] <- list()
genes2process <- unique(esscore_tab_outlier_drug$GENE[esscore_tab_outlier_drug$cancer == cancer_tmp & esscore_tab_outlier_drug$regulated])
genes2process
pairs2plot[[cancer_tmp]][["MET"]] <- c("MET:MET:Y1234") # is a known drug for kidney cancer
pairs2plot[[cancer_tmp]][["BRAF"]] <- c("BRAF:BRAF:S446") # S446 activating
pairs2plot[[cancer_tmp]][["MAP2K1"]] <- c("MAP2K1:MAPK1:Y187") # MAP2K1 have downstream effects and have been studied for RCC
pairs2plot[[cancer_tmp]][["MAPK1"]] <- c("MAPK1:STMN1:S25") # kidney cancer specific kinase-substrate pair

# bussiness ------------------------------------------------------------------
# for (cancer in c("BRCA")) {
for (cancer in cancers2process) {
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
    partIDs <- colnames(pho_tab)[!(colnames(pho_tab) %in% c("Gene", "Phosphosite", "Peptide_ID"))]
    sample_annotation <- data.frame(partID = partIDs)
    
    rna_tab <- loadRNA(cancer = cancer)
    maf <- loadMaf(cancer = cancer, maf_files = maf_files)
    if (any(maf$Hugo_Symbol %in% names(pairs2plot[[cancer]]))) {
      mut_mat <- generate_somatic_mutation_matrix(pair_tab = names(pairs2plot[[cancer]]), maf = maf)
    } else {
      mut_mat <- maf %>%
        filter(Hugo_Symbol %in% names(pairs2plot[[cancer]]))
    }
    
    mut_mat4sample_anno <- generate_somatic_mutation_matrix(pair_tab = mut_genes2w, maf = maf)
    mut_mat4sample_anno <- t(mut_mat4sample_anno)
    mut_df4sample_anno <- data.frame(mut_mat4sample_anno)
    mut_df4sample_anno$partID <- rownames(mut_df4sample_anno)
    sample_annotation <- merge(sample_annotation, mut_df4sample_anno, by = c("partID"), all.x = T)
    
    cna_tab <- loadCNAstatus(cancer = cancer)
    
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
    
    
    # add subtypes to column annotation ---------------------------------------
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
      df_value <- NULL
      
      # input Mutation ----------------------------------------------------------
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
          mut_mat2merge <- data.frame(partID = partIDs, expression_type = "Mutation", row_id = paste0(geneA, "_Mutation"))
          mut_mat2merge$color_value <- ifelse(mut_mat2merge$partID %in% mut_mat.m$partID[mut_mat.m$Gene == gene & mut_mat.m$variant_class_sim != "wild_type"], expressionType2colorValue(vector_expressionType = "Mutation"), 0)
          df_value <- rbind(df_value, mut_mat2merge)
        }
      } else {
        print("no mutation!")
      }
      
      
      # input CNA data ----------------------------------------------------------
      cna_tab_tmp <- cna_tab[cna_tab$gene %in% geneA, ]
      if (nrow(cna_tab_tmp) > 0) {
        cna_tab.m <- melt(cna_tab_tmp, id.vars = "gene")
        colnames(cna_tab.m) <- c("gene", "partID", "CNA")
        
        for (gene in unique(cna_tab.m$gene[cna_tab.m$CNA != "neutral"])) {
          cna_mat2merge <- data.frame(partID = partIDs, expression_type = "CNA", row_id = paste0(geneA, "_CNA"))
          cna_mat2merge$color_value <- ifelse(cna_mat2merge$partID %in% cna_tab.m$partID[cna_tab.m$gene == gene & cna_tab.m$CNA == "amplification"], expressionType2colorValue(vector_expressionType = "CNA"), 0)
          df_value <- rbind(df_value, cna_mat2merge)
        }
      } else {
        print("no CNA!")
      }
      
      # input RNA data ------------------------------------------------------
      rna_tab_tmp <- rna_tab %>%
        filter(gene == geneA)
      if (!is.null(rna_tab_tmp)) {
        if (nrow(rna_tab_tmp) > 0) {
          rna_tab_scaled <- rna_tab_tmp[, intersect(colnames(rna_tab_tmp), partIDs)]
          rna_tab_scaled <- scale_by_row(rna_tab_scaled)
          rna_outlier_tmp <- (rna_tab_scaled > outlier_sd)
          rna_outlier.m <- melt(rna_outlier_tmp)
          rna_outlier.m$color_value <- rna_outlier.m$value
          rna_outlier.m$color_value[!is.na(rna_outlier.m$value)] <- ifelse(rna_outlier.m$color_value[!is.na(rna_outlier.m$value)] == TRUE, expressionType2colorValue(vector_expressionType = "RNA"), 0)
          rna_outlier2merge <- rna_outlier.m %>%
            mutate(row_id = paste0(geneA, "_RNA")) %>%
            mutate(partID = Var2) %>%
            mutate(expression_type = "RNA") %>%
            select(partID, expression_type, row_id, color_value)
          df_value <- rbind(df_value, rna_outlier2merge)
        }
      }
      
      # input protein data ------------------------------------------------------
      pro_outlier_tmp <- escore_tab_outlier[escore_tab_outlier$SELF == "cis" & escore_tab_outlier$GENE == geneA,]
      pro_outlier_tmp <- pro_outlier_tmp[1,]
      if (!is.null(pro_outlier_tmp)) {
        if (nrow(pro_outlier_tmp) > 0) {
          pro_outlier.m <- melt(pro_outlier_tmp, 
                                id.vars = colnames(pro_outlier_tmp)[!(colnames(pro_outlier_tmp) %in% partIDs)])
          pro_outlier.m$color_value <- pro_outlier.m$value
          pro_outlier.m$color_value[!is.na(pro_outlier.m$value)] <- ifelse(pro_outlier.m$color_value[!is.na(pro_outlier.m$value)] == TRUE, expressionType2colorValue(vector_expressionType = "PRO"), 0)
          pro_outlier2merge <- pro_outlier.m %>%
            mutate(row_id = paste0(geneA, "_PRO")) %>%
            mutate(partID = variable) %>%
            mutate(expression_type = "PRO") %>%
            select(partID, expression_type, row_id, color_value)
          df_value <- rbind(df_value, pro_outlier2merge)
        }
      }
      
      # Input phosphoprotein data -----------------------------------------------
      phog_outlier_tmp <- escore_tab_outlier[escore_tab_outlier$SELF == "trans" & escore_tab_outlier$GENE == geneA,]
      phog_outlier_tmp <- phog_outlier_tmp[1,]
      if (!is.null(phog_outlier_tmp)) {
        if (nrow(phog_outlier_tmp) > 0) {
          phog_outlier.m <- melt(phog_outlier_tmp, 
                                 id.vars = colnames(phog_outlier_tmp)[!(colnames(phog_outlier_tmp) %in% partIDs)])
          phog_outlier.m$color_value <- phog_outlier.m$value
          phog_outlier.m$color_value[!is.na(phog_outlier.m$value)] <- ifelse(phog_outlier.m$color_value[!is.na(phog_outlier.m$value)] == TRUE, expressionType2colorValue(vector_expressionType = "PHO"), 0)
          phog_outlier2merge <- phog_outlier.m %>%
            mutate(row_id = paste0(geneA, "_PHO")) %>%
            mutate(partID = variable) %>%
            mutate(expression_type = "PHO") %>%
            select(partID, expression_type, row_id, color_value)
          
          df_value <- rbind(df_value, phog_outlier2merge)
        }
      }
      
      # input kinase-substrate score outliers ----------------------------
      esscore_tab_outlier_tmp <- esscore_tab_outlier[esscore_tab_outlier$pair %in% pairs,]
      if (!is.null(esscore_tab_outlier_tmp)) {
        if (nrow(esscore_tab_outlier_tmp) > 0) {
          esscore_tab_outlier.m <- melt(esscore_tab_outlier_tmp[, c("pair", "SUB_GENE", "SUB_MOD_RSD", partIDs)], 
                                        id.vars = c("pair", "SUB_GENE", "SUB_MOD_RSD"))
          esscore_tab_outlier.m <- esscore_tab_outlier.m %>%
            mutate(site_id = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
            select(variable, value, pair, site_id)
          
          site_ids_ordered <- sort(unique(esscore_tab_outlier.m$site_id))
          site_ids_ordered <- c(site_ids_ordered[!grepl(x = site_ids_ordered, pattern = geneA)], site_ids_ordered[grepl(x = site_ids_ordered, pattern = geneA)])
          esscore_tab_outlier.m$site_id <- factor(esscore_tab_outlier.m$site_id, levels = site_ids_ordered)
          esscore_tab_outlier.m <- esscore_tab_outlier.m[order(esscore_tab_outlier.m$site_id),]
          
          if (nrow(esscore_tab_outlier.m) > 0) {
            for (pair_tmp in unique(esscore_tab_outlier.m$pair[!is.na(esscore_tab_outlier.m$value) & esscore_tab_outlier.m$value == T])) {
              esscore_tab_outlier2merge <- esscore_tab_outlier.m %>%
                filter(pair == pair_tmp)
              esscore_tab_outlier2merge$color_value <- esscore_tab_outlier2merge$value
              esscore_tab_outlier2merge$color_value[!is.na(esscore_tab_outlier2merge$value)] <- ifelse(esscore_tab_outlier2merge$color_value[!is.na(esscore_tab_outlier2merge$value)] == TRUE, expressionType2colorValue(vector_expressionType = "KS"), 0)
              
              esscore_tab_outlier2merge <- esscore_tab_outlier2merge %>%
                mutate(row_id = paste0(geneA, "_KS_", paste0(str_split(string = pair, pattern = ":")[[1]], collapse = "_"))) %>%
                mutate(partID = variable) %>%
                mutate(expression_type = "KS") %>%
                select(partID, expression_type, row_id, color_value)
              df_value <- rbind(df_value, esscore_tab_outlier2merge)
            }
          }
        }
      }
      
      df_value <- unique(df_value)
      df_value <- df_value %>%
        mutate(id = paste0(row_id, ":", partID)) %>%
        filter(!duplicated(id)) %>%
        unique
      
      # make the matrix for the heatmap body ------------------------------------
      mat_value <- dcast(data = df_value, row_id ~ partID, value.var = "color_value")
      mat_value %>% head()
      rownames(mat_value) <- mat_value$row_id
      mat_value <- as.matrix(mat_value[,-1])
      
      mat_values_sup <- rbind(mat_values_sup, mat_value)
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
    
    
    my_heatmap <- pheatmap(mat_values_sup, 
                           color = c("#d9d9d9", rainbow(length(unique(df_value$color_value[!is.na(df_value$color_value) & df_value$color_value > 0])) + 2)),
                           gaps_row = gaps_row_vector,
                           na_col = "white", cellwidth = 4, cellheight = 5,
                           cluster_rows=F, cluster_cols=F, fontsize_row = 6,
                           number_color = "black", fontsize_number = 15,
                           show_colnames = F, annotation_colors = ann_colors)
    
    save_pheatmap_pdf(x = my_heatmap, 
                      filename = fn, 
                      width = 10, height = 20)
    
    
    # add expression outlier status -------------------------------------------
    expression_outlier_mat4sample_anno <- t(mat_values_sup)
    expression_outlier_mat4sample_anno[!is.na(expression_outlier_mat4sample_anno) & expression_outlier_mat4sample_anno != 0] <- TRUE
    expression_outlier_mat4sample_anno <- data.frame(expression_outlier_mat4sample_anno)
    expression_outlier_mat4sample_anno$partID <- rownames(expression_outlier_mat4sample_anno)
    sample_annotation <- merge(sample_annotation, expression_outlier_mat4sample_anno, by = c("partID"), all.x = T)
    
    # input sample mapping file -----------------------------------------------
    sample_map <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/A6_PGDAC_Data_Clear_Cell_Renal_Cell_Carcinoma/MountSinai/CPTAC3_CCRCC/Batch1-4/CPTAC-3 CCRCC samples tumor-normal info.csv", data.table = F)
    sample_map <- unique(sample_map[, c("Aliquot ID", "Case ID", "Tissue Type")])
    rownames(sample_map) <- sample_map$`Aliquot ID`
    sample_map <- sample_map[!(sample_map$`Case ID` %in% c("C3N-01175", "C3N-01180", "C3N-00313", "C3N-00435", "C3N-00832")),]
    
    # input ESTIMATE result ---------------------------------------------------
    estimate_tab <- readxl::read_excel("./Ding_Lab/Projects_Current/CPTAC/CPTACIII/CCRCC_AWG/CCRCC_shared_data/manuscripts/Table S7.xlsx", sheet = "ESTIMATE scores")
    estimate_purity_rna <- as.numeric(as.vector(estimate_tab[7,2:176]))
    estimate_purity_rna_df <- data.frame(ESTIMATE_TumorPurity_RNA = estimate_purity_rna, sampID = unlist(estimate_tab[3,2:176]))
    estimate_purity_rna_df$partID <- sample_map[estimate_purity_rna_df$sampID, "Case ID"]
    estimate_purity_rna_df$tissue_type <- sample_map[estimate_purity_rna_df$sampID, "Tissue Type"]
    estimate_purity_rna_df <- estimate_purity_rna_df %>%
      filter(tissue_type == "Tumor")
    estimate_purity_pro <- as.numeric(as.vector(estimate_tab[18,2:184]))
    estimate_purity_pro_df <- data.frame(ESTIMATE_TumorPurity_PRO = estimate_purity_pro, sampID = unlist(estimate_tab[14,2:184]))
    estimate_purity_pro_df$partID <- sample_map[estimate_purity_pro_df$sampID, "Case ID"]
    estimate_purity_pro_df$tissue_type <- sample_map[estimate_purity_pro_df$sampID, "Tissue Type"]
    estimate_purity_pro_df <- estimate_purity_pro_df %>%
      filter(tissue_type == "Tumor")
    sample_annotation <- merge(sample_annotation, estimate_purity_rna_df %>%
                                 select(partID, ESTIMATE_TumorPurity_RNA), by = c("partID"), all.x = T)
    sample_annotation <- merge(sample_annotation, estimate_purity_pro_df %>%
                                 select(partID, ESTIMATE_TumorPurity_PRO), by = c("partID"), all.x = T)
    sample_annotation <- merge(sample_annotation, sample_map[sample_map$`Tissue Type` == "Tumor",], by.x = c("partID"), by.y = c("Case ID"), all.x = T)
    
    # input xcell result and add immnue group ------------------------------------------------------
    xcell_tab <- readxl::read_excel("./Ding_Lab/Projects_Current/CPTAC/CPTACIII/CCRCC_AWG/CCRCC_shared_data/manuscripts/Table S7.xlsx", sheet = "xCell Signatures", skip = 2)
    ## immnue group 
    immune_groups <- unlist(xcell_tab[1,2:ncol(xcell_tab)])
    immune_group_df <- data.frame(immune_group = immune_groups, sampID = names(immune_groups))
    sample_annotation <- merge(sample_annotation, immune_group_df, by.x = c("Aliquot ID"), by.y = c("sampID"), all.x = T)
    sample_annotation <- sample_annotation %>%
      filter(!is.na(immune_group))
    
    
    # add xcell estimate of cell types ----------------------------------------
    xcell_tab$Samples
    xcell_cell_types <- c("")
    
    # add cyber sort estiamte of cell types -----------------------------------
    
    
    # add arm level CNV from GISTIC2 ------------------------------------------
    
    
    
    
    
    sample_annotation_fn <- paste0(subdir1, "CPTAC3_ccRCC_discovery_set_sample_annotation.csv")
    write.table(x = sample_annotation, file = sample_annotation_fn, sep = ",", row.names = F, col.names = T, quote = F)
  }
}
