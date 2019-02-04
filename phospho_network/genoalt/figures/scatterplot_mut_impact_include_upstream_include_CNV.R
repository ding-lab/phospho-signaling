# Yige Wu @ WashU 2018 Apr
## check regulated pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/expression_matrices/expression_shared.R')

library(ggrepel)
library(readxl)
library(ggpubr)

# set variables -----------------------------------------------------------
my_comparisons <- list( enzyme_mutation = c("enzyme_mutation", "control"))
amp_thres <- log2(1.3)
del_thres <- log2(0.7)

# plot mutation impact -----------------------------------------------------------
resultDnow <- makeOutDir(resultD = resultD)

## decide on the kinase and substrate
enzyme <- "AKT1"; substrate <- "GSK3B"; rsd <- "S9"; enzyme_upstream <- c("ERBB2", "EGFR", "IGF1R", "PTEN", "PIK3R1", "PAK1", "PIK3CA")
enzyme <- "AKT1"; substrate <- "BAD"; rsd <- "S99"; enzyme_upstream <- c("ERBB2", "EGFR", "IGF1R", "PTEN", "PIK3R1", "PAK1", "PIK3CA")
enzyme <- "PIK3CA"; substrate <- "AKT1"; rsd <- "T308"; enzyme_upstream <- c("ERBB2", "EGFR", "IGF1R", "PTEN", "PIK3R1")
enzyme <- "AKT1"; substrate <- "BAD"; rsd <- "S134"; enzyme_upstream <- c("ERBB2", "EGFR", "IGF1R", "PTEN", "PIK3R1", "PAK1", "PIK3CA")

## make directory
subdir1 <- paste0(resultDnow, enzyme, "/")
dir.create(subdir1)
subdir2 <- paste0(subdir1, substrate, "/")
dir.create(subdir2)
subdir3 <- paste0(subdir2, rsd, "/")
dir.create(subdir3)

## do for each cancer type
for (cancer in c("BRCA", "CO")) {
  subdir4 <- paste0(subdir3, cancer, "/")
  dir.create(subdir4)
  
  ## get the expression level of the substrate
  if (rsd == "PRO") {
    affected_exp_data <- loadProteinNormalizedTumor(cancer = cancer)
    colnames(affected_exp_data)[1] <- "Gene"
    affected_exp_head <- data.frame(SUBSTRATE = affected_exp_data$Gene)
    affected_exp_substrate <- affected_exp_data[affected_exp_head$SUBSTRATE == substrate,]
  } else {
    affected_exp_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
    affected_exp_head <- data.frame(str_split_fixed(string = affected_exp_data$Phosphosite, pattern = ":", n = 2))
    colnames(affected_exp_head) <- c("SUBSTRATE", "SUB_MOD_RSD")
    affected_exp_head$SUB_MOD_RSD <- toupper(affected_exp_head$SUB_MOD_RSD)
    affected_exp_head$SUBSTRATE <- affected_exp_data$Gene
    affected_exp_substrate <- affected_exp_data[affected_exp_head$SUBSTRATE == substrate & affected_exp_head$SUB_MOD_RSD == rsd,]
  }
  ## get sample ID map
  samples <- colnames(affected_exp_data)[!(colnames(affected_exp_data) %in% c("Gene", "Phosphosite")) & !(colnames(affected_exp_data) %in% c("idx")) & !grepl(pattern = "Mix", x = (colnames(affected_exp_data)))  & !grepl(pattern = "rep", x = (colnames(affected_exp_data)))]
  sampIDs <- samples
  sample_map <- loadSampMap()
  partIDs <- sampID2partID(sampleID_vector = samples, sample_map = sample_map)
  names(partIDs) <- samples
  affected_exp_substrate <- affected_exp_substrate[, sampIDs]
  colnames(affected_exp_substrate) <- partIDs
  
  ## get mutation matrix
  mut_mat <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_impact_proteome/", cancer, "_somatic_mutation.txt"), data.table = F)
  
  ## get CNV table
  cnv <- fread(input = paste0(cptac2p_genomicD, "copy_number/gatk/v1.3.swap_contamination_fixed/prospective_somatic/gene_level/v1.3.CPTAC2_prospective.2018-03-19/", toupper(substr(cancer, start = 1, stop = 2)), "/gene_level_CNV.", substr(cancer, start = 1, stop = 2), ".v1.3.2018-03-19.tsv"), data.table = F)
  
  partIDs_overlap <- intersect(partIDs, colnames(mut_mat))
  partIDs_overlap <- intersect(colnames(cnv), intersect(partIDs, colnames(mut_mat)))
  length(partIDs_overlap)
  
  ## get the patient IDs with enzyme alterations
  enzyme_mut_partIDs <- partIDs_overlap[(mut_mat[mut_mat$Hugo_Symbol == enzyme, partIDs_overlap] != "") & (mut_mat[mut_mat$Hugo_Symbol == enzyme, partIDs_overlap] != "Silent")]
  enzyme_cnv_partIDs <- partIDs_overlap[(cnv[cnv$gene == enzyme, partIDs_overlap] > amp_thres) | (cnv[cnv$gene == enzyme, partIDs_overlap] < del_thres)]
  
  ## get the patient IDs with enzyme upstream alterations
  upstream_mut_partIDs <- partIDs_overlap[colSums((mut_mat[mut_mat$Hugo_Symbol %in% enzyme_upstream, partIDs_overlap] != "") & (mut_mat[mut_mat$Hugo_Symbol %in% enzyme_upstream, partIDs_overlap] != "Silent")) > 0]
  upstream_cnv_partIDs <- partIDs_overlap[colSums((cnv[cnv$gene %in% enzyme_upstream, partIDs_overlap] > amp_thres) | (cnv[cnv$gene %in% enzyme_upstream, partIDs_overlap] < del_thres)) > 0]
  
  ## get the patient IDs with substrate alterations
  if (substrate %in% mut_mat$Hugo_Symbol) {
    substrate_mut_partIDs <- partIDs_overlap[(mut_mat[mut_mat$Hugo_Symbol == substrate, partIDs_overlap] != "") & (mut_mat[mut_mat$Hugo_Symbol == substrate, partIDs_overlap] != "Silent")]
  } else {
    substrate_mut_partIDs <- NULL
  }
  
  ## get the patient IDs of controls
  control_partIDs <- partIDs_overlap[!(partIDs_overlap %in% c(enzyme_mut_partIDs, upstream_mut_partIDs, substrate_mut_partIDs))]
  
  ## reshape and group by who had alterations
  df <- melt(affected_exp_substrate)
  colnames(df) <- c("partID", "sub_exp")
  df <- df[!is.na(df$sub_exp),]
  
  df$group <- "control"
  df$group[df$partID %in% c(upstream_mut_partIDs, upstream_cnv_partIDs)] <- "upstream_mutation"
  df$group[df$partID %in% c(enzyme_mut_partIDs, enzyme_cnv_partIDs)] <- "enzyme_mutation"
  df$group[df$partID %in% substrate_mut_partIDs] <- "substrate_mutation"
  
  ## boxplot
  ### get the 95% CI of controls
  substrate_exp_CI_low <- quantile(df$sub_exp[df$group == "control"], probs = 0.05, na.rm = T)
  substrate_exp_CI_high <- quantile(df$sub_exp[df$group == "control"], probs = 0.95, na.rm = T)
  
  ### mark samples with mutations who are outside the CI
  df$is.outlier_high <- (df$sub_exp > substrate_exp_CI_high)
  df$is.outlier_low <- (df$sub_exp < substrate_exp_CI_low)
  
  ## take a look at outlier mutations
  outlier_high_mut_mat <- mut_mat[mut_mat$Hugo_Symbol %in% c(enzyme, substrate, enzyme_upstream), df$partID[!is.na(df$is.outlier_high) & df$is.outlier_high]]
  outlier_low_mut_mat <- mut_mat[mut_mat$Hugo_Symbol %in% c(enzyme, substrate, enzyme_upstream), df$partID[!is.na(df$is.outlier_low) & df$is.outlier_low]]
  
  ## add text
  maf <- loadMaf(cancer = cancer, maf_files = maf_files)
  maf$partID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  maf$text <- paste0(maf$Hugo_Symbol, "_", maf$HGVSp_Short)
  df_text <- vector(mode = "character", length = nrow(df))
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  for (partID in df$partID[!is.na(df$is.outlier_high)]) {
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
  tab2p$x <- factor(tab2p$x, levels = c("upstream_mutation", "enzyme_mutation", "control"))
  tab2p$y <- as.vector(tab2p$sub_exp)
  
  
  pos <- position_jitter(width = 0.5, seed = 1)
  p = ggplot(tab2p, aes(x=x, y=y))
  p = p + geom_point(aes(color = group, shape = group), position = pos, stroke = 0, alpha = 0.6, size = 4)
  p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
  p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = group, label= text, color = group),
                          force = 1, segment.size = 0.5, segment.alpha = 0.2, size = 4, alpha=0.6, position = pos)
  p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio)"))
  p = p + theme_nogrid()
  p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 15, face = "bold"),
                axis.text.x = element_text(size= 10, vjust=0.5, hjust = 0.5, face = "bold"),
                axis.text.y = element_text(colour="black", size=8))
  p = p + theme(title = element_text(size = 18, face = "bold"))
  p = p + ggtitle(label = paste0(enzyme, " mutational association on ", substrate, "_", rsd))
  p
  fn = paste0(subdir4, enzyme, "_", substrate, "_", rsd, ".pdf")
  ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
  
  fn = paste0(subdir4, enzyme, "_", substrate, "_", rsd, ".png")
  ggsave(file=fn, height=7, width = 8, device = png())
  dev.off()
  
  ## do statisitical test
  
}





