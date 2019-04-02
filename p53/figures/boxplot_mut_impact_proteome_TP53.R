# Yige Wu @ WashU 2018 Apr
## check regulated pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source("./cptac2p_analysis/p53/TP53_shared.R")


library(readxl)
# set variables -----------------------------------------------------------
my_comparisons <- list(enzyme_mutation = c("enzyme_mutation", "control"))
sample_type2test <- names(my_comparisons)
cohort <- "TCGA"
# cohort <- "CPTAC"

# plot mutation impact -----------------------------------------------------------
resultDnow <- makeOutDir(resultD = resultD)
# enzyme <- "TP53"; substrate <- "ESR1"; rsd <- "RNA"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "ESR1"; rsd <- "PRO"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "MSH2"; rsd <- "PRO"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "MSH6"; rsd <- "PRO"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "MSH6"; rsd <- "RNA"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "MSH2"; rsd <- "RNA"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "BCL2"; rsd <- "PRO"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "AKT1"; rsd <- "PRO"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "BAK1"; rsd <- "PRO"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "BAK1"; rsd <- "RNA"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "STAT1"; rsd <- "RNA"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "STAT1"; rsd <- "PRO"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "PRKAB1"; rsd <- "RNA"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "PRKAB1"; rsd <- "PRO"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "PRKAB1"; rsd <- "PRO"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "CCNB1"; rsd <- "PRO"; enzyme_upstream <-  NULL
enzyme <- "TP53"; substrate <- "CCNB1"; rsd <- "RNA"; enzyme_upstream <-  NULL

## make directory
subdir1 <- paste0(resultDnow, enzyme, "/")
dir.create(subdir1)
subdir2 <- paste0(subdir1, substrate, "/")
dir.create(subdir2)
subdir3 <- paste0(subdir2, rsd, "/")
dir.create(subdir3)

## do for each cancer type
# for (cancer in c("BRCA", "UCEC", "LIHC", "CO", "CCRCC", "OV")) {
for (cancer in c("BRCA", "UCEC", "CO")) {
  if (cohort == "CPTAC") {
    ## get the expression level of the substrate
    if (cancer %in% c("BRCA", "OV", "CO")) {
      pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
      pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
    } else if (cancer == "UCEC") {
      pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
      pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
    } else if (cancer == "CCRCC") {
      pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
      pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
    } else if (cancer == "LIHC"){
      pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
      pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
    }
    rna_tab <- loadRNA(cancer = cancer)
    maf <- loadMaf(cancer = cancer, maf_files = maf_files)
    cna_tab <- loadCNAstatus(cancer = cancer)
  }
  if (cohort == "TCGA") {
    pro_tab <- loadTCGARPPA(cancer = cancer, expression_type = "PRO")
    pho_tab <- loadTCGARPPA(cancer = cancer, expression_type = "PHO")
    rna_tab <- loadTCGARNA(cancer = cancer)
    maf <- loadTCGAMaf(cancer = cancer)
    cna_tab <- loadTCGACNA(cancer = cancer)
  }
  
  if (rsd == "PRO") {
    affected_exp_data <- pro_tab
    affected_exp_head <- data.frame(SUBSTRATE = affected_exp_data$Gene)
    affected_exp_substrate <- affected_exp_data[affected_exp_head$SUBSTRATE == substrate & !is.na(affected_exp_head$SUBSTRATE),]
    affected_exp_substrate$Gene <- NULL
  } else if (rsd == "RNA") {
    affected_exp_data <- rna_tab
    affected_exp_head <- data.frame(SUBSTRATE = affected_exp_data$gene)
    affected_exp_substrate <- affected_exp_data[affected_exp_head$SUBSTRATE == substrate,]
    affected_exp_substrate$gene <- NULL
  } else {
    affected_exp_data <- pho_tab
    affected_exp_head <- data.frame(SUBSTRATE = affected_exp_data$Gene, SUB_MOD_RSD = affected_exp_data$Phosphosite)
    affected_exp_substrate <- affected_exp_data[affected_exp_head$SUBSTRATE == substrate & affected_exp_head$SUB_MOD_RSD == rsd  & !is.na(affected_exp_head$SUBSTRATE),]
    affected_exp_substrate$Gene <- NULL
    affected_exp_substrate$Phosphosite <- NULL
    affected_exp_substrate$Peptide_ID <- NULL
  }
  
  if (nrow(affected_exp_substrate) == 0 ){
    next()
  }
  
  partIDs <- colnames(affected_exp_data)[!(colnames(affected_exp_data) %in% c("Gene", "gene", "Phosphosite", "Peptide_ID"))]
  
  ## get mutation matrix
  pair_tab <- data.frame(GENE = enzyme, SUB_GENE = substrate)
  mut_mat <- generate_somatic_mutation_matrix(pair_tab = pair_tab, maf = maf)
  partIDs_in_mut <- colnames(mut_mat)[!(colnames(mut_mat) %in% "Hugo_Symbol")]
  if (cohort == "TCGA") {
    partID_tmp <- str_split_fixed(partIDs_in_mut, pattern = "-", n = 4)
    partIDs_in_mut <- paste(partID_tmp[,1], partID_tmp[,2], partID_tmp[,3], sep = "-")
  }
  colnames(mut_mat) <- c("Hugo_Symbol", partIDs_in_mut)
  
  # get the patient IDs with enzyme upstream mutations ----------------------
  if (!is.null(enzyme_upstream)) {
    upstream_mut_partIDs <- partIDs_in_mut[colSums((mut_mat[mut_mat$Hugo_Symbol %in% enzyme_upstream, partIDs_in_mut] != "") & (mut_mat[mut_mat$Hugo_Symbol %in% enzyme_upstream, partIDs_in_mut] != "Silent")) > 0]
    if (length(upstream_mut_partIDs) < 4 ){
      next()
    }
  } else {
    upstream_mut_partIDs <- NULL
  }
  
  ## get the patient IDs with substrate alterations
  if (substrate %in% mut_mat$Hugo_Symbol) {
    substrate_mut_partIDs <- partIDs_in_mut[(mut_mat[mut_mat$Hugo_Symbol == substrate, partIDs_in_mut] != "") & (mut_mat[mut_mat$Hugo_Symbol == substrate, partIDs_in_mut] != "Silent")]
  } else {
    substrate_mut_partIDs <- NULL
  }
  other_mut_partIDs <- partIDs_in_mut[!grepl(x = mut_mat[, partIDs_in_mut], pattern = "Missense_Mutation") & !grepl(x = mut_mat[, partIDs_in_mut], pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del") & (!is.na(mut_mat[, partIDs_in_mut]) &  mut_mat[, partIDs_in_mut] != "")]
  other_mut_partIDs <- other_mut_partIDs[!is.na(other_mut_partIDs)]
  missense_partIDs <- partIDs_in_mut[grepl(x = mut_mat[, partIDs_in_mut], pattern = "Missense_Mutation")]
  truncation_partIDs <- partIDs_in_mut[grepl(x = mut_mat[, partIDs_in_mut], pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del")]
  
  # get the patient IDs with enzyme upstream CNAs ----------------------
  cna_tab <- cna_tab[cna_tab$gene == enzyme,]
  cna_tab <- cna_tab[1,]
  partIDs_in_cna <- colnames(cna_tab)[!(colnames(cna_tab) %in% c("gene"))]
  if (cohort == "TCGA") {
    partID_tmp <- str_split_fixed(partIDs_in_cna, pattern = "-", n = 4)
    partIDs_in_cna <- paste(partID_tmp[,1], partID_tmp[,2], partID_tmp[,3], sep = "-")
  }
  colnames(cna_tab) <- c("gene", partIDs_in_cna)
  
  deletion_partIDs <- partIDs_in_cna[cna_tab[, partIDs_in_cna] == "deletion"]
  
  
  ## get the patient IDs of controls
  control_partIDs <- partIDs[!(partIDs %in% c(missense_partIDs, truncation_partIDs, deletion_partIDs,
                                              upstream_mut_partIDs, substrate_mut_partIDs))]
  control_partIDs
  
  ## reshape and group by who had alterations
  sup_tab <- melt(affected_exp_substrate)
  head(sup_tab)
  if (cohort == "CPTAC") {
    colnames(sup_tab) <- c("partID", "sub_exp")
  }
  if (cohort == "TCGA") {
    colnames(sup_tab) <- c("partID_exp", "sub_exp")
    partID_tmp <- str_split_fixed(sup_tab$partID_exp, pattern = "-", n = 4)
    sup_tab$partID <- paste(partID_tmp[,1], partID_tmp[,2], partID_tmp[,3], sep = "-")
    
  }
  
  
  # group samples by mutation status ----------------------------------------
  sup_tab$group <- "WT"
  sup_tab$group[sup_tab$partID %in% other_mut_partIDs] <- "Other"
  sup_tab$group[sup_tab$partID %in% deletion_partIDs] <- "Deletion"
  sup_tab$group[sup_tab$partID %in% truncation_partIDs] <- "Truncation"
  sup_tab$group[sup_tab$partID %in% missense_partIDs] <- "Missense"
  sup_tab$group[sup_tab$partID %in% truncation_partIDs & sup_tab$partID %in% deletion_partIDs] <- "Deletion_Truncation"
  sup_tab$group[sup_tab$partID %in% missense_partIDs & sup_tab$partID %in% deletion_partIDs] <- "Deletion_Missense"
  sup_tab$group[sup_tab$partID %in% missense_partIDs & sup_tab$partID %in% truncation_partIDs] <- "Missense_Truncation"
  
  stop("meow")
  
  table(sup_tab$group)
  # if (enzyme != substrate) {
  #   sup_tab$group[sup_tab$partID %in% substrate_mut_partIDs] <- "substrate_mutation"
  # }
  
  ## boxplot
  ### get the 95% CI of controls
  substrate_exp_CI_low <- quantile(sup_tab$sub_exp[sup_tab$group == "control"], probs = 0.1, na.rm = T)
  substrate_exp_CI_high <- quantile(sup_tab$sub_exp[sup_tab$group == "control"], probs = 0.9, na.rm = T)
  
  ### mark samples with mutations who are outside the CI
  sup_tab$is.outlier_high <- (sup_tab$sub_exp > substrate_exp_CI_high)
  sup_tab$is.outlier_low <- (sup_tab$sub_exp < substrate_exp_CI_low)
  
  ## take a look at outlier mutations
  outlier_high_mut_mat <- mut_mat[mut_mat$Hugo_Symbol %in% c(enzyme, substrate, enzyme_upstream), intersect(partIDs_in_mut, sup_tab$partID[!is.na(sup_tab$is.outlier_high) & sup_tab$is.outlier_high])]
  outlier_low_mut_mat <- mut_mat[mut_mat$Hugo_Symbol %in% c(enzyme, substrate, enzyme_upstream), intersect(partIDs_in_mut, sup_tab$partID[!is.na(sup_tab$is.outlier_low) & sup_tab$is.outlier_low])]
  
  ## add text
  maf$partID <- maf$Tumor_Sample_Barcode
  maf$text <- paste0(maf$Hugo_Symbol, "_", maf$HGVSp_Short)
  sup_tab_text <- vector(mode = "character", length = nrow(sup_tab))
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  for (partID in partIDs_in_mut) {
    mut_mat_tmp <- mut_mat[mut_mat$Hugo_Symbol %in% c(enzyme, substrate, enzyme_upstream), c("Hugo_Symbol", partID)]
    mut_genes <- rownames(mut_mat_tmp)[!(mut_mat_tmp[, partID] %in% c("", "Silent"))]
    if (length(mut_genes) > 0) {
      ## get the all the mutations for these genes for this patient
      maf_tmp <- maf$text[maf$partID == partID & (maf$Hugo_Symbol %in% mut_genes) & maf$Variant_Classification != "Silent"]
      maf_tmp <- maf$text[(maf$partID == partID) & (maf$Hugo_Symbol %in% mut_genes)]
      
      text_tmp <- paste0(maf_tmp, collapse = "\n")
      sup_tab_text[sup_tab$partID == partID] <- text_tmp
    }
  }
  sup_tab$text <- sup_tab_text
  
  ### make plot
  tab2p <- sup_tab
  tab2p$x <- tab2p$group
  tab2p$x <- sort_mutation_class_complex(tab2p$group)
  tab2p$y <- as.vector(tab2p$sub_exp)
  ## filtering
  tab2p <- tab2p[tab2p$x %in% names(colors_mut),]
  tab2p$y <- remove_outliers_IQR(x = tab2p$sub_exp, out_thres = 3, na.rm = T)
  tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
  
  
  subdir4 <- paste0(subdir3, cancer, "/")
  dir.create(subdir4)
  fn = paste0(subdir4, cohort, "_", cancer, "_", enzyme, "_", substrate, "_", rsd, "_with_upstream", paste0(enzyme_upstream, collapse = "_"), ".pdf")
  
  pos <- position_jitter(width = 0.2, seed = 1)
  
  if (cohort == "CPTAC") {
    p = ggplot(tab2p, aes(x=x, y=y))
    p = p + geom_violin(aes(fill = x),  color = NA, alpha = 1)
    p = p + geom_boxplot(width=.1)
    p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.6, size = 1)
    p = p + scale_fill_manual(values = colors_mut, name = "TP53 mutation status")
    if (rsd == "PRO") {
      p = p + labs(y=paste0(substrate, " protein abundance\n(log2 ratio)"))
    }
    if (rsd != "PRO" & rsd != "RNA") {
      p = p + labs(y=paste0(substrate, " ", rsd, "\nphosphorylation abundance (log2 ratio)"))
    }
    if (rsd == "RNA") {
      p = p + labs(y=paste0(substrate, " mRNA abundance"))
    }
    p = p + theme_bw() + theme_nogrid()
    p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 15, face = "bold"),
                  axis.text.x = element_blank(),
                  axis.text.y = element_text(colour="black", size=8), legend.title = element_text(size = 10))
    p = p + theme(title = element_text(size = 18, face = "bold"))
    p = p + stat_compare_means(data = tab2p, mapping = aes(x = x, y = y, label = ..p.signif..), symnum.args = symnum.args, ref.group = "WT", hide.ns = T)     # Add global Anova p-value
    p = p + scale_x_discrete(breaks = c("Deletion_Truncation", "Deletion_Missense", "Truncation", "Missense", "Deletion", "WT"))
    p
    ggsave(file=fn, height= 3, width = 6.5, useDingbats=FALSE)
    
  } else {
    p = ggplot(tab2p, aes(x=x, y=y))
    p = p + geom_violin(aes(fill = x),  color = NA, alpha = 1)
    p = p + geom_boxplot(width=.1)
    p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.6, size = 1)
    p = p + scale_fill_manual(values = colors_mut, name = "TP53 mutation status")
    if (rsd == "PRO") {
      p = p + labs(y=paste0(substrate, " protein abundance\n(RPPA)"))
    }
    if (rsd != "PRO" & rsd != "RNA") {
      p = p + labs(y=paste0(substrate, " ", rsd, "\nphosphorylation abundance (RPPA)"))
    }
    if (rsd == "RNA") {
      p = p + labs(y=paste0(substrate, " mRNA abundance"))
    }
    p = p + theme_bw() + theme_nogrid()
    p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 15, face = "bold"),
                  axis.text.x = element_blank(),
                  axis.text.y = element_text(colour="black", size=8), legend.title = element_text(size = 10))
    p = p + theme(title = element_text(size = 18, face = "bold"))
    p = p + stat_compare_means(data = tab2p, mapping = aes(x = x, y = y, label = ..p.signif..), symnum.args = symnum.args, ref.group = "WT", hide.ns = T)     # Add global Anova p-value
    p = p + scale_x_discrete(breaks = c("Deletion_Truncation", "Deletion_Missense", "Truncation", "Missense", "Deletion", "WT"))
    p
    ggsave(file=fn, height= 3, width = 6.5, useDingbats=FALSE)
  }
  
  # ## test enzyme_mutation to control
  # if (length(tab2p$sub_exp[tab2p$group == "enzyme_mutation" & !is.na(tab2p$sub_exp)]) == 0) {
  #   next()
  # }
  # print(wilcox.test(tab2p$sub_exp[tab2p$group == "enzyme_mutation" & !is.na(tab2p$sub_exp)], tab2p$sub_exp[tab2p$group == "control"& !is.na(tab2p$sub_exp)]))
  # 
  # ## do statisitical test on enrichment
  # if (enzyme == "CTNNB1") {
  #   tmp <- str_split_fixed(string = tab2p$text, pattern = "_p.[A-Z]", 2)[,2]
  #   tmp <- str_split_fixed(string = tmp, pattern = "[A-Z]|[a-z]|\\*", n = 2)[,1]
  #   tmp
  #   tab2p$position <- as.numeric(tmp)
  #   tab2p$is.GSK3B_target <- (tab2p$position %in% c(33,37,41))
  #   tab2p$is.CK1_target <- (tab2p$position %in% c(45))
  #   tab2p$is.GSK3B_region <- (tab2p$position <= 42 & tab2p$position >= 32)
  #   tab2p$is.GSK3B_CK1_region <- (tab2p$position <= 46 & tab2p$position >= 32)
  #   
  #   head(tab2p)
  #   tab2test_GSK3B_target <- table(tab2p[!is.na(tab2p$position) & !is.na(tab2p$sub_exp) & tab2p$group == "enzyme_mutation", c("is.GSK3B_target", "is.outlier_high")])
  #   fisher.test(tab2test_GSK3B_target)
  #   
  #   tab2test_GSK3B_region <- table(tab2p[!is.na(tab2p$position) & !is.na(tab2p$sub_exp) & tab2p$group == "enzyme_mutation", c("is.GSK3B_region", "is.outlier_high")])
  #   tab2test_GSK3B_region
  #   fisher.test(tab2test_GSK3B_region)
  #   
  #   tab2test_CK1_target <- table(tab2p[!is.na(tab2p$position) & !is.na(tab2p$sub_exp) & tab2p$group == "enzyme_mutation", c("is.CK1_target", "is.outlier_high")])
  #   fisher.test(tab2test_CK1_target)
  #   
  #   tab2test_GSK3B_CK1_region <- table(tab2p[!is.na(tab2p$position) & !is.na(tab2p$sub_exp) & tab2p$group == "enzyme_mutation", c("is.GSK3B_CK1_region", "is.outlier_high")])
  #   tab2test_GSK3B_CK1_region
  #   fisher.test(tab2test_GSK3B_region)
  #   
  #   table(tab2p[tab2p$group == "enzyme_mutation", c("is.outlier_high", "is.GSK3B_CK1_region")])
  # }
  # table(tab2p[tab2p$group == "enzyme_mutation", c("is.outlier_high")])
  
}





