# Yige Wu @ WashU 2018 Jan
# plotting bubble chart and subtype phosphorylation/protein expression heatmap sid_rowe by sid_rowe for BRCA dataset
## TODO: keep kinase on one side and keep substrate on the other side
## filter down the trans pairs

# library -----------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
library(grid_row)
library(dplyr)
library(grid_rowExtra)
library(gtable)
library(readxl)
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')

source(paste0(baseD, "cptac2p_analysis/pan3can_aes.R")) # aes for general purposes; it should be one directory out of the working directory

# functions ---------------------------------------------------------------
plotpheatmap <- function(mat_value, color.palette, col_anno, row_anno, ann_colors, fn, width, height) {
  result <- tryCatch({
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno, annotation_row = row_anno,
                           # scale = "row",
                           na_col = "white", 
                           gaps_row = (1:length(unique(row_anno$expression_type)))*length(unique(row_anno$subtype)),
                           cluster_cols=T, show_colnames = T, 
                           cluster_rows=F, show_rownames = T, annotation_colors = ann_colors)
  }, error = function(err) {
    print(print(paste("MY_ERROR:  ",err)))
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno, annotation_row = row_anno,
                           # scale = "row",
                           na_col = "white", 
                           gaps_row = (1:length(unique(row_anno$expression_type)))*length(unique(row_anno$subtype)),
                           cluster_cols=F, show_colnames = T, 
                           cluster_rows=F, show_rownames = T, annotation_colors = ann_colors)    
    return(NA)
  }, warning = function(war) {
    print(print(paste("MY_WARNING:  ", war)))
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno, annotation_row = row_anno,
                           # scale = "row",
                           na_col = "white", 
                           gaps_row = (1:length(unique(row_anno$expression_type)))*length(unique(row_anno$subtype)),
                           cluster_cols=F, show_colnames = T, 
                           cluster_rows=F, show_rownames = T, annotation_colors = ann_colors)    
    return(NA)
    
    return(NULL)
  }, finally = {
    save_pheatmap_pdf(x = my_heatmap, 
                      filename = fn, 
                      width = width, height = height)
  })
  return(result)
}

partID2subtype <- function(partIDs, cancer) {
  if (cancer == "BRCA") {
    subtypes <- partID2pam50(patientID_vector = partIDs, pam50_map = loadPAM50Map())
  } else if (cancer == "CO") {
    subtypes <- partID2MSI(patientID_vector = partIDs, subtype_map = loadMSIMap())
  } else if (cancer == "UCEC") {
    subtypes <- partID2UCECsubtype(patientID_vector = partIDs)
  } else {
    next()
  }
  return(subtypes)
}

# set variables -------------------------------------------------------------------
resultDnow <- makeOutDir(resultD = resultD)
var_list <- c('pro_kin', 'pho_kin')
names(var_list) <- c('cis', 'trans')
subtypes_list <- c('pam50', 'MSI_type', 'tumor_normal_WBasal'); names(subtypes_list) <- c("BRCA", "CO", "OV")

## plotting paramters
cap <- 3
breaks = seq(-(cap),cap, by=0.2)
## add color palette
color.palette <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(length(breaks))
color.palette.Rd <- col_paletteR(length(breaks))
## CHANGE THIS
# cancer <- "BRCA"
data2process <- matrix(data = c("BRCA", "CDAP", "tumor", "scaled", "cptac2p",
                                # "CO", "CDAP", "tumor", "scaled", "cptac2p",
                                "UCEC", "PGDAC", "tumor", "median_polishing", "cptac3"), ncol = 5, byrow = T)

# input table to filter for pairs to show ---------------------------------
fisher_stat_cancers <- fread(input = paste0(ppnD,"es_score/table/calculate_es_pair_score/fisher_stat_tab_outlier1.5SD_all_cancers_kinase_reg_nonNA20.txt"), data.table = F)


# Gene list to filter kinase-substrate pairs ------------------------------
genes2plot <- unique(c(unlist(SMGs), "EGFR", "ERBB2", "PKN1", "CDK1", "PAK2"))

# input regression result to show correlation -----------------------------
## load regression table
reg_nonNA <- 20
enzyme_type <- "kinase"
regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                   "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)



# BRCA PAM50 subtypes -----------------------------------------------------
for (i in 1:2) {
  cancer <- data2process[i,1]
  pipeline_type <- data2process[i,2]
  sample_type <- data2process[i,3]
  norm_type <- data2process[i,4]
  cptac_phase <- data2process[i,5]
  
  # input subtype aggregated protein/phosphoprotein/phosphosite data -------------
  # pho_subtype_mean <- fread(input = paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_",subtypes_list[cancer],"_mean.txt"), data.table = F)
  # pro_subtype_mean <- fread(input = paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_",subtypes_list[cancer],"_mean.txt"), data.table = F)
  # phog_subtype_mean <- fread(input = paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_collapsed_",subtypes_list[cancer],"_mean.txt"), data.table = F)
  # pho_subtype_mean$SUB_MOD_RSD <- toupper(as.vector(str_split_fixed(string = pho_subtype_mean$Phosphosite, pattern = ":", n = 2)[,2]))
  ## input gene-level and site-level phosphorylation and protein levels
  pro_data <- loadParseProteomicsData(cancer = cancer, expression_type = "PRO", sample_type = sample_type, pipeline_type = pipeline_type, norm_type = norm_type)
  
  phog_data <- loadParseProteomicsData(cancer = cancer, expression_type = "collapsed_PHO", sample_type = sample_type, pipeline_type = pipeline_type, norm_type = norm_type)
  
  pho_data <- loadParseProteomicsData(cancer = cancer, expression_type = "PHO", sample_type = sample_type, pipeline_type = pipeline_type, norm_type = norm_type)
  pho_data <- pho_data[, colnames(pho_data)[!(colnames(pho_data) %in% "Peptide_ID")]]
  
  for (SELF in c("cis", "trans")) {
    ## filter the enrichment results
    if (SELF == "cis") {
      fisher_stat_cancers2plot <- fisher_stat_cancers[fisher_stat_cancers$cancer == cancer & fisher_stat_cancers$p.value < 0.05 & (fisher_stat_cancers$GENE %in% genes2plot | fisher_stat_cancers$SUB_GENE %in% genes2plot) & fisher_stat_cancers$SELF == SELF,]
      fisher_stat_cancers2plot <- fisher_stat_cancers2plot[!(fisher_stat_cancers2plot$GENE == "ERBB2" & !(fisher_stat_cancers2plot$pair %in% c("ERBB2:ERBB2:Y1109", "ERBB2:ERBB2:Y1218"))),]
      fisher_stat_cancers2plot <- fisher_stat_cancers2plot[order(fisher_stat_cancers2plot$p.value),]
      fisher_stat_cancers2plot <- fisher_stat_cancers2plot[order(fisher_stat_cancers2plot$p.value),]
      fisher_stat_cancers2plot <- fisher_stat_cancers2plot[order(fisher_stat_cancers2plot$SELF),]
      pairs2plot <- unique(fisher_stat_cancers2plot$pair)
    } else {
      fisher_stat_cancers2plot <- fisher_stat_cancers[fisher_stat_cancers$cancer == cancer & fisher_stat_cancers$p.value < 0.05 & fisher_stat_cancers$SELF == SELF,]
      fisher_stat_cancers2plot2add <- fisher_stat_cancers[(fisher_stat_cancers$GENE == "CDK1" & fisher_stat_cancers$SUB_GENE %in% c("TP53", "TP53BP1", "DNMT1", "PTTG", "BUB1B", "MCM4", "LIG1", "LIG3", "RFC1", "ANAPC4", "RB1", "CUX1")),]
      fisher_stat_cancers2plot2add <- rbind(fisher_stat_cancers2plot2add, fisher_stat_cancers[(fisher_stat_cancers$GENE == "DYRK1A" & fisher_stat_cancers$SUB_GENE %in% c("CCND1")),])
      fisher_stat_cancers2plot2add <- rbind(fisher_stat_cancers2plot2add, fisher_stat_cancers[(fisher_stat_cancers$GENE == "MAPK9" & fisher_stat_cancers$SUB_GENE %in% c("JUN")),])
      fisher_stat_cancers2plot <- fisher_stat_cancers2plot2add
    }
    pairs2plot <- unique(fisher_stat_cancers2plot$pair)
    
    # construct the data matrix part of heatmap -------------------------------
    sup_tab <- NULL
    ## the heatmap columns will be each pair
    for (pair in pairs2plot) {
      geneA <- str_split(string = pair, pattern = ":")[[1]][1]
      geneB <- str_split(string = pair, pattern = ":")[[1]][2]
      phosphosite <- str_split(string = pair, pattern = ":")[[1]][3]
      
      ## add protein data
      if (geneA %in% pro_data$Gene) {
        pro_geneA <- pro_data[pro_data$Gene == geneA,]
        pro_geneA.m <- melt(pro_geneA, id.vars = c("Gene"))
        pro_geneA.m$subtype <- partID2subtype(partIDs = pro_geneA.m$variable, cancer = cancer)
        ## take mean averaged across subtypes and then scale to -1 to 1 (?)
        pro_geneA.m_subtype <- data.frame(subtype = unique(pro_geneA.m$subtype))
        pro_geneA.m_subtype$exp_value <- sapply(X = pro_geneA.m_subtype$subtype, FUN = function(subtype, mat) mean(mat$value[mat$subtype == subtype], na.rm = T), mat = pro_geneA.m)
        pro2merge <- data.frame(pair = pair, subtype = pro_geneA.m_subtype$subtype, exp_value = pro_geneA.m_subtype$exp_value, id_row = paste0("kinase_protein_", pro_geneA.m_subtype$subtype))
        sup_tab <- rbind(sup_tab, pro2merge)
      } 
      
      ## add protein data
      if (geneB %in% pro_data$Gene & SELF == "trans") {
        pro_geneB <- pro_data[pro_data$Gene == geneB,]
        pro_geneB.m <- melt(pro_geneB, id.vars = c("Gene"))
        pro_geneB.m$subtype <- partID2subtype(partIDs = pro_geneB.m$variable, cancer = cancer)
        ## take mean averaged across subtypes and then scale to -1 to 1 (?)
        pro_geneB.m_subtype <- data.frame(subtype = unique(pro_geneB.m$subtype))
        pro_geneB.m_subtype$exp_value <- sapply(X = pro_geneB.m_subtype$subtype, FUN = function(subtype, mat) mean(mat$value[mat$subtype == subtype], na.rm = T), mat = pro_geneB.m)
        pro2merge <- data.frame(pair = pair, subtype = pro_geneB.m_subtype$subtype, exp_value = pro_geneB.m_subtype$exp_value, id_row = paste0("substrate_protein_", pro_geneB.m_subtype$subtype))
        sup_tab <- rbind(sup_tab, pro2merge)
      } 
      
      ## add phosphoprotein data
      if (geneA %in% phog_data$Gene & SELF == "trans") {
        phog_geneA <- phog_data[phog_data$Gene == geneA,]
        phog_geneA.m <- melt(phog_geneA, id.vars = c("Gene"))
        phog_geneA.m$subtype <- partID2subtype(partIDs = phog_geneA.m$variable, cancer = cancer)
        ## take mean averaged across subtypes and then scale to -1 to 1 (?)
        phog_geneA.m_subtype <- data.frame(subtype = unique(phog_geneA.m$subtype))
        phog_geneA.m_subtype$exp_value <- sapply(X = phog_geneA.m_subtype$subtype, FUN = function(subtype, mat) mean(mat$value[mat$subtype == subtype], na.rm = T), mat = phog_geneA.m)
        phog2merge <- data.frame(pair = pair, subtype = phog_geneA.m_subtype$subtype, exp_value = phog_geneA.m_subtype$exp_value, id_row = paste0("kinase_phosphoprotein_", phog_geneA.m_subtype$subtype))
        sup_tab <- rbind(sup_tab, phog2merge)
      } 
      
      ## add phosphosite data
      pho_geneB <- pho_data[pho_data$Gene == geneB & pho_data$Phosphosite == phosphosite,]
      if (nrow(pho_geneB) > 0) {
        pho_geneB.m <- melt(pho_geneB, id.vars = c("Gene", "Phosphosite"))
        pho_geneB.m$subtype <- partID2subtype(partIDs = pho_geneB.m$variable, cancer = cancer)
        ## take mean averaged across subtypes and then scale to -1 to 1 (?)
        pho_geneB.m_subtype <- data.frame(subtype = unique(pho_geneB.m$subtype))
        pho_geneB.m_subtype$exp_value <- sapply(X = pho_geneB.m_subtype$subtype, FUN = function(subtype, mat) mean(mat$value[mat$subtype == subtype], na.rm = T), mat = pho_geneB.m)
        pho2merge <- data.frame(pair = pair, subtype = pho_geneB.m_subtype$subtype, exp_value = pho_geneB.m_subtype$exp_value, id_row = paste0("substrate_phosphosite_", pho_geneB.m_subtype$subtype))
        sup_tab <- rbind(sup_tab, pho2merge)
      }
    }
    sup_tab <- sup_tab[sup_tab$subtype != "no_RNA-seq",]
    sup_tab <- sup_tab[order(sup_tab$subtype, sup_tab$id_row),]
    
    ## make the matrix for the heatmap body
    df_value <- dcast(data = sup_tab, id_row ~ pair, value.var = "exp_value")
    df_value %>% head()
    mat_value <- as.matrix(df_value[,-1])
    rownames(mat_value) <- df_value$id_row
    head(mat_value)
    mat_value %>% head()
    
    # construct the column annotation part ------------------------------------
    col_anno <- data.frame(pair = pairs2plot)
    ann_colors <- list()
    
    row_anno <- data.frame(id_row = unique(sup_tab$id_row))
    tmp <- str_split_fixed(string = row_anno$id_row, pattern = "_", n = 3)
    row_anno$expression_type <- paste0(tmp[,1], "_", tmp[,2])
    row_anno$subtype <- tmp[,3]
    rownames(row_anno) <- row_anno$id_row
    row_anno$id_row <- NULL
    
    ## add correlation coeffcient from regression results
    regression2merge <- regression[regression$pair %in% pairs2plot & regression$Cancer == cancer,]
    regression2merge$corr_coef_beta <- regression2merge$coef_pro_kin
    regression2merge$corr_coef_beta[regression2merge$SELF == "trans"] <- regression2merge$coef_pho_kin[regression2merge$SELF == "trans"]
    
    if (nrow(regression2merge) > 0) {
      col_anno <- merge(col_anno, regression2merge[, c("pair", "corr_coef_beta")], all.x = T)
      tmp <- color.palette.Rd; 
      names(tmp) <- seq(from = -1, to = 1, length.out = length(color.palette.Rd))
      ann_colors[["corr_coef_beta"]] <- tmp
    }
    
    # actual plotting ---------------------------------------------------------
    rownames(col_anno) <- col_anno$pair
    col_anno$pair <- NULL
    
    fn <- paste0(resultDnow, SELF, "_ks_pairs_in_", cancer, "subtypes.pdf")
    plotpheatmap(mat_value, color.palette, col_anno, row_anno, ann_colors, fn, ifelse(SELF == "cis", 6, 9), ifelse(SELF == "cis", 5, 5))
  }
  
}
