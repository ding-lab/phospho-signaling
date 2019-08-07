# Yige Wu @ WashU 2018 Apr
## check regulated pairs

# source ------------------------------------------------------------------
# source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
# source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_plotting.R")
source(path2phospho_network_shared)
source('./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/expression_matrices/expression_shared.R')




library(ggrepel)
library(readxl)
library(ggpubr)

# set variables -----------------------------------------------------------
my_comparisons <- list(enzyme_mutation = c("enzyme_mutation", "control"))
sample_type2test <- names(my_comparisons)


# plot mutation impact -----------------------------------------------------------
resultDnow <- makeOutDir()

## decide on the kinase and substrate
# enzyme <- "CTNNB1"; substrate <- "CDK6"; rsd <- "PRO"; enzyme_upstream <-  c("APC")
# enzyme <- "CTNNB1"; substrate <- "PAK1"; rsd <- "T212"; enzyme_upstream <-  c("APC")
# enzyme <- "CTNNB1"; substrate <- "CSNK1A1"; rsd <- "PRO"; enzyme_upstream <-  c("APC")
# enzyme <- "AKT1"; substrate <- "GSK3B"; rsd <- "S9"; enzyme_upstream <-  c("PIK3CA")
# enzyme <- "AKT1"; substrate <- "BAD"; rsd <- "S99"; enzyme_upstream <-  c("PIK3CA")
# enzyme <- "TP53"; substrate <- "CHEK2"; rsd <- "T432"; enzyme_upstream <-  NULL
# enzyme <- "TP53"; substrate <- "CHEK2"; rsd <- "PRO"; enzyme_upstream <-  NULL
# enzyme <- "CTNNB1"; substrate <- "PRKD1"; rsd <- "PRO"; enzyme_upstream <-  NULL
# enzyme <- "CTNNB1"; substrate <- "PRKD1"; rsd <- "S205"; enzyme_upstream <-  "APC"
# enzyme <- "APC"; substrate <- "GSK3B"; rsd <- "T390"; enzyme_upstream <-  c("CTNNB1", "FBXW7")
# enzyme <- "CTNNB1"; substrate <- "CTNNB1"; rsd <- "S191"; enzyme_upstream <-  "APC"
# enzyme <- "CTNNB1"; substrate <- "CTNNB1"; rsd <- "S552"; enzyme_upstream <-  NULL
# enzyme <- "CTNNB1"; substrate <- "CTNNB1"; rsd <- "PRO"; enzyme_upstream <-  "APC"
# enzyme <- "APC"; substrate <- "CTNNB1"; rsd <- "PRO"; enzyme_upstream <-  NULL
# enzyme <- "APC"; substrate <- "CTNNB1"; rsd <- "PRO"; enzyme_upstream <-  c("CTNNB1")
# enzyme <- "APC"; substrate <- "PRKD1"; rsd <- "S205"; enzyme_upstream <-  c("CTNNB1")
# enzyme <- "APC"; substrate <- "CTNNB1"; rsd <- "S191"; enzyme_upstream <-  c("CTNNB1")
# enzyme <- "FBXW7"; substrate <- "GSK3B"; rsd <- "T390"; enzyme_upstream <-  c("CTNNB1", "APC")
# enzyme <- "FBXW7"; substrate <- "CTNNB1"; rsd <- "S191"; enzyme_upstream <-  c("CTNNB1", "APC")
# enzyme <- "FBXW7"; substrate <- "CTNNB1"; rsd <- "PRO"; enzyme_upstream <-  c("CTNNB1", "APC")
# enzyme <- "TP53"; substrate <- "IKBKB"; rsd <- "PRO"; enzyme_upstream <- NULL
# enzyme <- "JAK1"; substrate <- "STAT3"; rsd <- "Y705"; enzyme_upstream <- NULL
# enzyme <- "TP53"; substrate <- "CDK1"; rsd <- "PRO"; enzyme_upstream <- NULL
# enzyme <- "TP53"; substrate <- "TP53BP1"; rsd <- "PRO"; enzyme_upstream <- NULL
# enzyme <- "TP53"; substrate <- "TP53BP1"; rsd <- "S1683"; enzyme_upstream <- NULL
# enzyme <- "TP53"; substrate <- "CDK1"; rsd <- "Y15"; enzyme_upstream <- NULL
# enzyme <- "TP53"; substrate <- "ATR"; rsd <- "PRO"; enzyme_upstream <- NULL
# enzyme <- "TP53"; substrate <- "ATR"; rsd <- "T1989"; enzyme_upstream <- NULL
# enzyme <- "MTOR"; substrate <- "EIF4EBP1"; rsd <- "S70"; enzyme_upstream <-  c("PIK3CA")
# enzyme <- "MTOR"; substrate <- "EIF4EBP1"; rsd <- "S65"; enzyme_upstream <-  c("PIK3CA")
# enzyme <- "AKT1"; substrate <- "LMNA"; rsd <- "S301"; enzyme_upstream <-  c("PIK3CA")
# enzyme <- "TP53"; substrate <- "CDK1"; rsd <- "RNA"; enzyme_upstream <- NULL
# enzyme <- "TP53"; substrate <- "ATR"; rsd <- "RNA"; enzyme_upstream <- NULL
# enzyme <- "KEAP1"; substrate <- "KEAP1"; rsd <- "PRO"; enzyme_upstream <- c("NFE2L2")
# enzyme <- "KEAP1"; substrate <- "NFE2L2"; rsd <- "PRO"; enzyme_upstream <- c("NFE2L2")
# enzyme <- "KEAP1"; substrate <- "NFE2L2"; rsd <- "S215"; enzyme_upstream <- c("NFE2L2")
# enzyme <- "KEAP1"; substrate <- "NFE2L2"; rsd <- "S433"; enzyme_upstream <- c("NFE2L2")
# pho_tab %>%
#   filter(Gene == "NFE2L2") %>%
#   select(Phosphosite) %>%
#   unique()
enzyme <- "PBRM1"; substrate <- "SMARCC1"; rsd <- "RNA"; enzyme_upstream <- NULL


  
## make directory
subdir1 <- paste0(resultDnow, enzyme, "/")
dir.create(subdir1)
subdir2 <- paste0(subdir1, substrate, "/")
dir.create(subdir2)
subdir3 <- paste0(subdir2, rsd, "/")
dir.create(subdir3)

## do for each cancer type
# for (cancer in c("BRCA", "UCEC", "LIHC", "CO", "OV", "CCRCC")) {
for (cancer in c( "CCRCC")) {
  subdir4 <- paste0(subdir3, cancer, "/")
  dir.create(subdir4)
  
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
  
  if (rsd == "PRO") {
    affected_exp_data <- pro_tab
    affected_exp_head <- data.frame(SUBSTRATE = affected_exp_data$Gene)
    affected_exp_substrate <- affected_exp_data[affected_exp_head$SUBSTRATE == substrate & !is.na(affected_exp_head$SUBSTRATE),]
    affected_exp_substrate$Gene <- NULL
  } else if (rsd == "RNA") {
    affected_exp_data <- loadRNA(cancer = cancer)
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

  partIDs <- colnames(affected_exp_data)
  ## get mutation matrix
  maf <- loadMaf(cancer = cancer, maf_files = maf_files)
  mut_mat <- generate_somatic_mutation_matrix(pair_tab = unique(c(enzyme, substrate, enzyme_upstream)), maf = maf)
  partIDs_overlap <- intersect(partIDs, colnames(mut_mat))
  partIDs_overlap <- partIDs
  
  ## get the patient IDs with enzyme alterations
  enzyme_mut_partIDs <- colnames(mut_mat)[(mut_mat[mut_mat$Hugo_Symbol == enzyme, ] != "") & (mut_mat[mut_mat$Hugo_Symbol == enzyme, ] != "Silent")]
  
  
  ## get the patient IDs with enzyme upstream alterations
  if (!is.null(enzyme_upstream)) {
    upstream_mut_partIDs <- colnames(mut_mat)[colSums((mut_mat[mut_mat$Hugo_Symbol %in% enzyme_upstream, ] != "") & (mut_mat[mut_mat$Hugo_Symbol %in% enzyme_upstream, ] != "Silent")) > 0]
    if (length(upstream_mut_partIDs) < 4 ){
      # next()
    }
  } else {
    upstream_mut_partIDs <- NULL
  }
  
  ## get the patient IDs with substrate alterations
  if (substrate %in% mut_mat$Hugo_Symbol) {
    substrate_mut_partIDs <- colnames(mut_mat)[(mut_mat[mut_mat$Hugo_Symbol == substrate, ] != "") & (mut_mat[mut_mat$Hugo_Symbol == substrate, ] != "Silent")]
  } else {
    substrate_mut_partIDs <- NULL
  }
  
  ## get the patient IDs of controls
  control_partIDs <- partIDs_overlap[!(partIDs_overlap %in% c(enzyme_mut_partIDs, upstream_mut_partIDs, substrate_mut_partIDs))]
  
  
  ## reshape and group by who had alterations
  df <- melt(affected_exp_substrate)
  head(df)
  colnames(df) <- c("partID", "sub_exp")
  if (is.null(enzyme_upstream)) {
    df <- df[df$partID %in% c(control_partIDs, enzyme_mut_partIDs),]
  }
  df$group <- "control"
  df$group[df$partID %in% upstream_mut_partIDs] <- "upstream_mutation"
  df$group[df$partID %in% enzyme_mut_partIDs] <- "enzyme_mutation"
  # if (enzyme != substrate) {
  #   df$group[df$partID %in% substrate_mut_partIDs] <- "substrate_mutation"
  # }
  
  ## boxplot
  ### get the 95% CI of controls
  substrate_exp_CI_low <- quantile(df$sub_exp[df$group == "control"], probs = 0.1, na.rm = T)
  substrate_exp_CI_high <- quantile(df$sub_exp[df$group == "control"], probs = 0.9, na.rm = T)
  
  ### mark samples with mutations who are outside the CI
  df$is.outlier_high <- (df$sub_exp > substrate_exp_CI_high)
  df$is.outlier_low <- (df$sub_exp < substrate_exp_CI_low)
  
  ## take a look at outlier mutations
  outlier_high_mut_mat <- mut_mat[mut_mat$Hugo_Symbol %in% c(enzyme, substrate, enzyme_upstream), intersect(colnames(mut_mat), df$partID[!is.na(df$is.outlier_high) & df$is.outlier_high])]
  outlier_low_mut_mat <- mut_mat[mut_mat$Hugo_Symbol %in% c(enzyme, substrate, enzyme_upstream), intersect(colnames(mut_mat), df$partID[!is.na(df$is.outlier_low) & df$is.outlier_low])]
  
  ## add text
  maf <- loadMaf(cancer = cancer, maf_files = maf_files)
  maf$partID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  maf$text <- paste0(maf$Hugo_Symbol, "_", maf$HGVSp_Short)
  df_text <- vector(mode = "character", length = nrow(df))
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  for (partID in colnames(mut_mat)) {
    mut_mat_tmp <- mut_mat[mut_mat$Hugo_Symbol %in% c(enzyme, substrate, enzyme_upstream), c("Hugo_Symbol", partID)]
    mut_genes <- rownames(mut_mat_tmp)[!(mut_mat_tmp[, partID] %in% c("", "Silent"))]
    if (length(mut_genes) > 0) {
      ## get the all the mutations for these genes for this patient
      maf_tmp <- maf$text[maf$partID == partID & (maf$Hugo_Symbol %in% mut_genes) & maf$Variant_Classification != "Silent"]
      maf_tmp <- maf$text[(maf$partID == partID) & (maf$Hugo_Symbol %in% mut_genes)]
      
      text_tmp <- paste0(maf_tmp, collapse = "\n")
      df_text[df$partID == partID] <- text_tmp
    }
  }
  df$text <- df_text
  
  ### make plot
  tab2p <- df
  tab2p$x <- tab2p$group
  tab2p$x <- factor(tab2p$x, levels = c( "enzyme_mutation", "upstream_mutation", "control",  "substrate_mutation"))
  tab2p$y <- as.vector(tab2p$sub_exp)
  tab2p <- unique(tab2p)
  pos <- position_jitter(width = 0.2, seed = 1)
  p = ggplot(tab2p, aes(x=x, y=y))
  p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.6, size = 2)
  p = p + geom_boxplot(aes(fill = group),  color = NA, alpha = 0.4, notch = T)
  p = p + scale_fill_manual(values = c("enzyme_mutation" = set1[1], "control" = "grey50", "upstream_mutation" = set1[1]))
  # p = p + geom_text_repel(data = tab2p[!is.na(df$is.outlier_high) & (df$is.outlier_high | df$is.outlier_low),], mapping = aes(segment.color = group, label= text, color = group),
  #                         force = 1, segment.size = 0.5, segment.alpha = 0.2, size = 2, alpha=0.6, position = pos)
  if (rsd == "PRO") {
    p = p + labs(y=paste0(substrate, " protein abundance(log2 ratio)"))
  } else if (rsd == "RNA") {
    p = p + labs(y=paste0(substrate, " RNA abundance(log2)"))
  } else {
    p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation\n(log2 ratio)"))
    
  }
  p = p + theme_nogrid()
  p = p + theme(title = element_text(size = 18, face = "bold"))
  p = p + stat_compare_means(data = tab2p, mapping = aes(x = x, y = y, label = ..p.signif..), symnum.args = symnum.args, ref.group = "control")     # Add global Anova p-value
  p = p + scale_x_discrete(breaks = c("upstream_mutation", 
                                      "enzyme_mutation", "substrate_mutation", "control"),
                           label = c(paste0(enzyme_upstream, "\nmutated"),
                                     paste0(enzyme, "\nmutated"),
                                     paste0(substrate, "\nmutated"),
                                     "control"))
  p <- p + guides(fill = F)
  if (!is.null(enzyme_upstream)) {
    p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10, face = "bold"),
                  axis.text.x = element_text(size= 15, vjust=0.5, hjust = 0.5, face = "bold"),
                  axis.text.y = element_text(colour="black", size=8))
    p
    fn = paste0(subdir4, enzyme, "_", substrate, "_", rsd, "_with_upstream", paste0(enzyme_upstream, collapse = "_"), ".pdf")
    ggsave(file=fn, height= 3, width = 4, useDingbats=FALSE)
  } else {
    p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10, face = "bold"),
                  axis.text.x = element_text(size= 20, vjust=0.5, hjust = 0.5, face = "bold"),
                  axis.text.y = element_text(colour="black", size=8))
    p
    fn = paste0(subdir4, enzyme, "_", substrate, "_", rsd, "_with_upstream", paste0(enzyme_upstream, collapse = "_"), ".pdf")
    ggsave(file=fn, height= 3, width = 3*length(unique(tab2p$group))/2, useDingbats=FALSE)
  }

# stop("")
  
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





