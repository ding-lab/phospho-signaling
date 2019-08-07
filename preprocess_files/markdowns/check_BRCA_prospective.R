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
regression %>% nrow()
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)
regression %>% nrow()

regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression <- annotate_ks_source(regression = regression)

# input outlier results -----------------------------------------------
esscore_tab_outlier_drug <- fread(input = "./cptac2p/analysis_results/phospho_network/druggability/figures/grid_percent_patient_with_druggable_outlier/esscore_tab_outlier_drug_genes.txt", data.table = F)
esscore_tab_outlier_drug <- annotate_ks_source(regression = esscore_tab_outlier_drug)
esscore_tab_outlier_drug <- esscore_tab_outlier_drug %>%
  mutate(pair_cancer = paste0(pair, ":", cancer)) %>%
  mutate(regulated = (pair_cancer %in% regression$pair_cancer[regression$regulated == T])) %>%
  mutate(regulated_uniq = (pair_cancer %in% regression$pair_cancer[regression$regulated_uniq == T]))

# inputs regulatory sites------------------------------------------------------------------
regulatory_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/Mar_04_2019/Regulatory_sites", data.table = F)
regulatory_sites <- regulatory_sites %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1])
regulatory_sites2merge <- regulatory_sites %>%
  select(GENE, SUB_MOD_RSD, ON_FUNCTION, ON_PROCESS, ON_PROT_INTERACT, ON_OTHER_INTERACT) %>%
  unique
esscore_tab_outlier_drug <- merge(esscore_tab_outlier_drug, regulatory_sites2merge, by.x = c("SUB_GENE", "SUB_MOD_RSD"), by.y = c("GENE", "SUB_MOD_RSD"), all.x = T)


pairs2plot <- list();
# input ccrcc -------------------------------------------------------------
cancer_tmp <- "CCRCC"
pairs2plot[[cancer_tmp]] <- list()
genes2process <- unique(esscore_tab_outlier_drug$GENE[esscore_tab_outlier_drug$cancer == cancer_tmp & esscore_tab_outlier_drug$regulated])
genes2process
pairs2plot[[cancer_tmp]][["MET"]] <- c("MET:MET:Y1234") # is a known drug for kidney cancer
pairs2plot[[cancer_tmp]][["BRAF"]] <- c("BRAF:BRAF:S446") # S446 activating
pairs2plot[[cancer_tmp]][["MAP2K1"]] <- c("MAP2K1:MAPK1:Y187") # MAP2K1 have downstream effects and have been studied for RCC
pairs2plot[[cancer_tmp]][["MAPK1"]] <- c("MAPK1:STMN1:S25") # kidney cancer specific kinase-substrate pair

# input CO -------------------------------------------------------------
# cancer_tmp <- "CO"
# esscore_tab_outlier_drug %>%
#   filter(cancer == cancer_tmp) %>%
#   filter(regulated_uniq == T)
# esscore_tab_outlier_drug %>%
#   filter(cancer == cancer_tmp) %>%
#   filter(regulated == T)
# pairs2plot[[cancer_tmp]] <- list()
# genes2process <- unique(esscore_tab_outlier_drug$GENE[esscore_tab_outlier_drug$cancer == cancer_tmp & esscore_tab_outlier_drug$regulated])
# genes2process
# pairs2plot[[cancer_tmp]][["PRKDC"]] <- c("PRKDC:XRCC6:S27", "PRKDC:PRKDC:S3205")
# #@ colorectal cancer specific, PRKDC:XRCC6:S27 is the only trans pair for PRKDC
# ## PRKDC:PRKDC:S3205 outlier overlap with PRKDC:XRCC6:S27 the most and have the highest coverage among the PRKDC cis pairs
# pairs2plot[[cancer_tmp]][["SRC"]] <- c("SRC:PTPRA:Y798")
# #@ colorectal cancer specific kinase-substrate pair, PTPRA:Y798 is direct target by SRC
# ## overlap with SRC copy number amplification, high data coverage
# ## PTPRA:Y798 induce PRPRA activity, implicated in cytoskeleton reorginization
# pairs2plot[[cancer_tmp]][["AKT1"]] <- c("AKT1:AKT1:S129") 
# ## having patients with both mutation and KS outlier score
# ## already needed to be measured for breast cancer samples
# ## cis-associated and have outlier KS score for breast, ovarian and colorectal cancer
# ## although for clinical trials only certain AKT inhibitor (GDC-0068) showed response in colorectal cancer for PTEN-deficient/PIK3CA-mutated colorectal cancer
# ## https://ascopubs.org/doi/10.1200/JCO.2011.35.5263
# ## A phase Ib study of the Akt inhibitor GDC-0068 with docetaxel (D) or mFOLFOX-6 (F) in patients (pts) with advanced solid tumors.
# ## A first-inhuman trial of GDC-0068: A novel, oral, ATP-competitive Akt inhibitor, demonstrates robust suppression of the Akt pathway in surrogate and tumor tissues.
# esscore_tab_outlier_drug %>%
#   filter(cancer == cancer_tmp) %>%
#   filter(regulated_uniq == T) %>%
#   filter(GENE == "SRC")
# 
# esscore_tab_outlier_drug %>%
#   filter(cancer == cancer_tmp) %>%
#   filter(SELF == "cis") %>%
#   filter(GENE == "SRC")
# 
# esscore_tab_outlier_drug %>%
#   filter(SELF == "cis") %>%
#   filter(pair == "AKT1:AKT1:S129")
# 
# regression %>%
#   filter(pair == "AKT1:AKT1:S129")
#   
# 
# cancers2process <- c("CO")

# input UCEC -------------------------------------------------------------
# cancer_tmp <- "UCEC"
# pairs2plot[[cancer_tmp]] <- list()
# genes2process <- unique(esscore_tab_outlier_drug$GENE[esscore_tab_outlier_drug$cancer == cancer_tmp & esscore_tab_outlier_drug$regulated])
# genes2process
# pairs2plot[[cancer_tmp]][["PLK1"]] <- c("PLK1:TP53BP1:S1623")
# ## PLK1 is an interesting drugable target for UCEC, TP53BP1 S1623 is direct target of PLK1
# ## from our mutational section we showed that TP53 mutations lead to high PLK1 and high TP53BP1 phospho
# ## TP53BP1 S1623 has been reported to inhibit DNA damage
# pairs2plot[[cancer_tmp]][["EGFR"]] <- c("EGFR:SHC1:Y428", "EGFR:EGFR:T693")
# ## EGFR is an interesting drugable target for UCEC, SHC1 Y428 is direct target of EGFR
# ## EGFR T693 has reported function of inducing interaction between EGF and EGFR, high coverage
# ## mutliple patients with overlapping outlier events in both pairs
# pairs2plot[[cancer_tmp]][["BRAF"]] <- c("BRAF:BRAF:S446")
# ## already going to validate in CCRCC
# ## overlaps with BRAF CNA and downstream MAPK1 activation, probably overlaps with KRAS mutations
# ## clinical trials looks promising for KRAS-mutated endometrial cancer: http://cancerres.aacrjournals.org/content/76/14_Supplement/CT005.short
# cancers2process <- c("UCEC")
# esscore_tab_outlier_drug %>%
#   filter(SELF == "cis") %>%
#   filter(pair == "AKT1:AKT1:S129")
# 
# esscore_tab_outlier_drug %>%
#   filter(cancer == cancer_tmp) %>%
#   filter(SELF == "cis") %>%
#   filter(GENE == "BRAF")

# input BRCA -------------------------------------------------------------
cancer_tmp <- "BRCA"
pairs2plot[[cancer_tmp]] <- list()
genes2process <- unique(esscore_tab_outlier_drug$GENE[esscore_tab_outlier_drug$cancer == cancer_tmp & esscore_tab_outlier_drug$regulated])
genes2process
pairs2plot[[cancer_tmp]][["AKT1"]] <- c("AKT1:GSK3B:S9", "AKT1:AKT1:S129")
## AKT1 is a SMG for BRCA
## AKT1 is an interesting drugable target for BRCA
## GSK3B S9 is direct target of AKT1, its phosphorylation inhibit apoptosis
## AKT1:GSK3B:S9 is breast cancer specific
## AKT1:S129 induce AKT1 activity, inhibit apoptosis and cis-associated
pairs2plot[[cancer_tmp]][["SRC"]] <- c("SRC:SRC:Y419", "SRC:PTPN6:Y536")
## SRC is an interesting drugable target for BRCA
## SRC:Y419 induce SRC activity
cancers2process <- c("BRCA")

esscore_tab_outlier_drug_brca <- esscore_tab_outlier_drug %>%
  filter(cancer == cancer_tmp) %>%
  filter(SELF == "trans") %>%
  filter(GENE == "SRC")

# input OV -------------------------------------------------------------
cancer_tmp <- "OV"
pairs2plot[[cancer_tmp]] <- list()
genes2process <- unique(esscore_tab_outlier_drug$GENE[esscore_tab_outlier_drug$cancer == cancer_tmp & esscore_tab_outlier_drug$regulated])
genes2process
pairs2plot[[cancer_tmp]][["SRC"]] <- c("SRC:CDK1:Y15")
## targeting SRC in preclinical studies: https://www.nature.com/articles/oncsis201654, https://www.nature.com/articles/cddis2017367
pairs2plot[[cancer_tmp]][["AKT1"]] <- c("AKT1:AKT1:S129")
## AKT inhibitor show durable tumor activity in ovarian cancer: https://academic.oup.com/annonc/article/28/suppl_5/mdx367.025/4108615
cancers2process <- c("OV")

esscore_tab_outlier_drug_ov <- esscore_tab_outlier_drug %>%
  filter(cancer == cancer_tmp) %>%
  filter(regulated == T)

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
    rna_tab <- loadRNA(cancer = cancer)
    maf <- loadMaf(cancer = cancer, maf_files = maf_files)

    if (any(maf$Hugo_Symbol %in% names(pairs2plot[[cancer]]))) {
      mut_mat <- generate_somatic_mutation_matrix(pair_tab = names(pairs2plot[[cancer]]), maf = maf)
    } else {
      mut_mat <- maf %>%
        filter(Hugo_Symbol %in% names(pairs2plot[[cancer]]))
    }
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
      
      # plot upset plot ---------------------------------------------------------
      tab2upset <- data.frame(partID = partIDs)
      for (expression_type_tmp in c("Mutation", "CNA", "RNA", "PRO", "PHO", "KS")) {
        tab2upset[, expression_type_tmp] <- as.numeric(tab2upset$partID %in% df_value$partID[(df_value$expression_type == expression_type_tmp & df_value$color_value > 0 & !is.na(df_value$color_value))])
      }
      
      sub_fn = paste(subdir1, geneA, "_", cancer , '_outlier.pdf',sep ="")
      pdf(sub_fn, height = 3, width = 4, onefile = F)
      upset(tab2upset, sets = intersect(c("Mutation", "CNA", "RNA", "PRO", "PHO", "KS"), unique(df_value$expression_type[df_value$color_value > 0 & !is.na(df_value$color_value)])), sets.bar.color = "#56B4E9", keep.order = T,
            empty.intersections = NULL)
      dev.off()
      

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
                           color = c("#d9d9d9", rainbow(length(unique(df_value$color_value[!is.na(df_value$color_value) & df_value$color_value > 0])) + 1)),
                           gaps_row = gaps_row_vector,
                           na_col = "white", cellwidth = 4, cellheight = 5,
                           cluster_rows=F, cluster_cols=F, fontsize_row = 6,
                           number_color = "black", fontsize_number = 15,
                           show_colnames = F, annotation_colors = ann_colors)
    
    save_pheatmap_pdf(x = my_heatmap, 
                      filename = fn, 
                      width = 10, height = 20)
    
  }
}

test <- regression %>%
  filter(Cancer == "BRCA") %>%
  filter(regulated_uniq == T) %>%
  filter(SUB_GENE %in% c("BRCA1", "BRCA2"))



## subtype info, need to adapt when subtype info is absent
# if ("subtype" %in% colnames(col_anno)) {
#   subtype_colors <- set2[1:length(unique(subtypes))]
#   names(subtype_colors) <- unique(subtypes)
#   ann_colors[["subtype"]] <- subtype_colors
# }

# if (("subtype" %in% colnames(col_anno))) {
#   col_anno <- merge(col_anno, subtypes2merge, by = c("partID"), all.x = T)
# } else {
#   col_anno$subtype <- ""
# }
# 
# rownames(col_anno) <- col_anno$partID
# col_anno$partID <- NULL
