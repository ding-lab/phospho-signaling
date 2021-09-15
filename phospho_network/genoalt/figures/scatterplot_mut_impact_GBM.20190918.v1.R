# Yige Wu @ WashU 2018 Apr
## for plotting trans mutatinal effect on proteome/phosphoproteome in GBM

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/phospho_network/phospho_network_plotting.R")

# set variables -----------------------------------------------------------
my_comparisons <- list(mut_gene_mutation = c("Missense_Group", "Control"))
sample_type2test <- names(my_comparisons)
cancer <- "GBM"

# input proteomics data ---------------------------------------------------
pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "d4")
pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "d6")
rna_tab <- fread(input = "./Ding_Lab/Projects_Current/CPTAC3-GBM/Data_Freeze/cptac3_gbm_data_freeze_v2.0.20190905/gene_expression/washu/rnaseq_washu_fpkm.v2.0.20190905.tsv.gz", data.table = F)
partIDs <- unique(c(colnames(pro_tab), colnames(pho_tab)))
rna_tab <- rna_tab[, c("gene_name", intersect(colnames(rna_tab), partIDs))]
colnames(rna_tab)[1] <- "gene"
rna_tab %>%
  head()
rna_mat <- as.matrix(rna_tab[,-1])
rna_mat_log2 <- log2(rna_mat+1)
rna_tab <- data.frame(gene = rna_tab$gene)
rna_tab <- cbind(rna_tab, as.data.frame(rna_mat_log2))
rna_tab %>%
  head()


# plot mutation impact -----------------------------------------------------------
mut_gene <- "TP53"; exp_gene <- "MDM2"; exp_site <- "PRO"; mut_gene_upstream <- NULL
mut_gene <- "TP53"; exp_gene <- "MDM2"; exp_site <- "RNA"; mut_gene_upstream <- NULL
mut_gene <- "TP53"; exp_gene <- "CDKN2A"; exp_site <- "RNA"; mut_gene_upstream <- NULL
mut_gene <- "TP53"; exp_gene <- "CDKN2A"; exp_site <- "PRO"; mut_gene_upstream <- NULL

## make directory
subdir1 <- paste0(makeOutDir(), mut_gene, "/")
dir.create(subdir1)
subdir2 <- paste0(subdir1, exp_gene, "/")
dir.create(subdir2)
subdir3 <- paste0(subdir2, exp_site, "/")
dir.create(subdir3)

if (exp_site == "PRO") {
  affected_exp_data <- pro_tab
  affected_exp_head <- data.frame(exp_gene = affected_exp_data$Gene)
  exp_df <- affected_exp_data[affected_exp_head$exp_gene == exp_gene & !is.na(affected_exp_head$exp_gene),]
  exp_df$Gene <- NULL
} else if (exp_site == "RNA") {
  affected_exp_data <- rna_tab
  affected_exp_head <- data.frame(exp_gene = affected_exp_data$gene)
  exp_df <- affected_exp_data[affected_exp_head$exp_gene == exp_gene,]
  exp_df$gene <- NULL
} else {
  affected_exp_data <- pho_tab
  affected_exp_head <- data.frame(exp_gene = affected_exp_data$Gene, SUB_MOD_exp_site = affected_exp_data$Phosphosite)
  exp_df <- affected_exp_data[affected_exp_head$exp_gene == exp_gene & affected_exp_head$SUB_MOD_exp_site == exp_site  & !is.na(affected_exp_head$exp_gene),]
  exp_df$Gene <- NULL
  exp_df$Phosphosite <- NULL
  exp_df$Peptide_ID <- NULL
}

if (nrow(exp_df) == 0 ){
  next()
}

partIDs <- colnames(affected_exp_data)
partIDs <- partIDs[!(partIDs %in% c("Gene", "Phosphosite"))]


## get mutation matrix
maf <- loadMaf(cancer = cancer, maf_files = maf_files)
mut_mat <- generate_somatic_mutation_matrix(pair_tab = unique(c(mut_gene, exp_gene, mut_gene_upstream)), maf = maf)

## get the patient IDs with mut_gene alterations
partIDs_overlap <- intersect(partIDs, colnames(mut_mat))
mut_partIDs <- partIDs_overlap[(mut_mat[mut_mat$Hugo_Symbol == mut_gene, partIDs_overlap] != "") & (mut_mat[mut_mat$Hugo_Symbol == mut_gene, partIDs_overlap] != "Silent")]
missense_partIDs <- partIDs_overlap[(mut_mat[mut_mat$Hugo_Symbol == mut_gene, partIDs_overlap] == "Missense_Mutation")]
truncation_partIDs <- partIDs_overlap[(mut_mat[mut_mat$Hugo_Symbol == mut_gene, partIDs_overlap] %in% c("Nonsense_Mutation", "Splice_Site", "In_Frame_Del", "Frame_Shift_Del", "Frame_Shift_Ins"))]
other_mut_partIDs <- mut_partIDs[!(mut_partIDs %in% c(missense_partIDs, truncation_partIDs))]
  
mut_partIDs

## get the patient IDs with mut_gene upstream alterations
if (!is.null(mut_gene_upstream)) {
  upstream_mut_partIDs <- partIDs_overlap[colSums((mut_mat[mut_mat$Hugo_Symbol %in% mut_gene_upstream, partIDs_overlap] != "") & (mut_mat[mut_mat$Hugo_Symbol %in% mut_gene_upstream, partIDs_overlap] != "Silent")) > 0]
  if (length(upstream_mut_partIDs) < 4 ){
    # next()
  }
} else {
  upstream_mut_partIDs <- NULL
}

## get the patient IDs with exp_gene alterations
if (exp_gene %in% mut_mat$Hugo_Symbol) {
  exp_gene_mut_partIDs <- partIDs_overlap[(mut_mat[mut_mat$Hugo_Symbol == exp_gene, partIDs_overlap] != "") & (mut_mat[mut_mat$Hugo_Symbol == exp_gene, partIDs_overlap] != "Silent")]
} else {
  exp_gene_mut_partIDs <- NULL
}

## get the patient IDs of controls
control_partIDs <- partIDs[!(partIDs %in% c(mut_partIDs, upstream_mut_partIDs, exp_gene_mut_partIDs))]
control_partIDs

## reshape and group by who had alterations
tab2p_sup <- melt(exp_df)

colnames(tab2p_sup) <- c("partID", "exp_value")
if (is.null(mut_gene_upstream)) {
  tab2p_sup <- tab2p_sup[tab2p_sup$partID %in% c(control_partIDs, mut_partIDs),]
}

# grouping samples for boxplot x axis --------------------------------------------------------
tab2p_sup$group <- "Control"
tab2p_sup$group[tab2p_sup$partID %in% upstream_mut_partIDs] <- "Upstream_Mut_Group"
tab2p_sup$group[tab2p_sup$partID %in% mut_partIDs] <- "Mut_Group"
tab2p_sup$group[tab2p_sup$partID %in% exp_gene_mut_partIDs] <- "Exp_Gene_Mut_Group"

tab2p_sup$group_mut_detailed <- "Control"
tab2p_sup$group_mut_detailed[tab2p_sup$partID %in% upstream_mut_partIDs] <- "Upstream_Mut_Group"
tab2p_sup$group_mut_detailed[tab2p_sup$partID %in% missense_partIDs] <- "Missense_Group"
tab2p_sup$group_mut_detailed[tab2p_sup$partID %in% truncation_partIDs] <- "Trunc_Group"
tab2p_sup$group_mut_detailed[tab2p_sup$partID %in% other_mut_partIDs] <- "Other_Mut_Group"
tab2p_sup$group[tab2p_sup$partID %in% exp_gene_mut_partIDs] <- "Exp_Gene_Mut_Group"

## boxplot
### get the 95% CI of controls
exp_gene_exp_CI_low <- quantile(tab2p_sup$exp_value[tab2p_sup$group == "control"], probs = 0.1, na.rm = T)
exp_gene_exp_CI_high <- quantile(tab2p_sup$exp_value[tab2p_sup$group == "control"], probs = 0.9, na.rm = T)

### mark samples with mutations who are outside the CI
tab2p_sup$is.outlier_high <- (tab2p_sup$exp_value > exp_gene_exp_CI_high)
tab2p_sup$is.outlier_low <- (tab2p_sup$exp_value < exp_gene_exp_CI_low)

## take a look at outlier mutations
outlier_high_mut_mat <- mut_mat[mut_mat$Hugo_Symbol %in% c(mut_gene, exp_gene, mut_gene_upstream), intersect(colnames(mut_mat), tab2p_sup$partID[!is.na(tab2p_sup$is.outlier_high) & tab2p_sup$is.outlier_high])]
outlier_low_mut_mat <- mut_mat[mut_mat$Hugo_Symbol %in% c(mut_gene, exp_gene, mut_gene_upstream), intersect(colnames(mut_mat), tab2p_sup$partID[!is.na(tab2p_sup$is.outlier_low) & tab2p_sup$is.outlier_low])]

## add mutation details
maf$partID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
maf$text <- paste0(maf$Hugo_Symbol, "_", maf$HGVSp_Short)
mut_text <- vector(mode = "character", length = nrow(tab2p_sup))
rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
for (partID in partIDs_overlap) {
  mut_mat_tmp <- mut_mat[mut_mat$Hugo_Symbol %in% c(mut_gene, exp_gene, mut_gene_upstream), c("Hugo_Symbol", partID)]
  mut_genes <- rownames(mut_mat_tmp)[!(mut_mat_tmp[, partID] %in% c("", "Silent"))]
  if (length(mut_genes) > 0) {
    ## get the all the mutations for these genes for this patient
    maf_tmp <- maf$text[maf$partID == partID & (maf$Hugo_Symbol %in% mut_genes) & maf$Variant_Classification != "Silent"]
    maf_tmp <- maf$text[(maf$partID == partID) & (maf$Hugo_Symbol %in% mut_genes)]
    
    text_tmp <- paste0(maf_tmp, collapse = "\n")
    mut_text[tab2p_sup$partID == partID] <- text_tmp
  }
}
tab2p_sup$text <- mut_text

# make plot for Detailed Mut vs WT -----------------------------------------
tab2p <- tab2p_sup
tab2p <- tab2p %>%
  filter(!is.na(exp_value))
tab2p$x <- tab2p$group_mut_detailed
tab2p$x <- factor(tab2p$x, levels = c( "Missense_Group", "Trunc_Group", "Other_Mut_Group",
                                       "Upstream_Mut_Group", "Exp_Gene_Mut_Group", "Control"))
tab2p$y <- as.vector(tab2p$exp_value)
tab2p <- unique(tab2p)
pos <- position_jitter(width = 0.2, seed = 1)
p = ggplot(tab2p, aes(x=x, y=y))
p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.6, size = 2)
p = p + geom_boxplot(aes(fill = group),  color = NA, alpha = 0.4, notch = T)
p = p + scale_fill_manual(values = c("Mut_Group" = set1[1], "Upstream_Mut_Group" = set1[1], "Exp_Gene_Mut_Group" = set1[1], "Control" = "grey50"))
if (exp_site == "PRO") {
  p = p + labs(y=paste0(exp_gene, " protein abundance(log2 ratio)"))
} else if (exp_site == "RNA") {
  p = p + labs(y=paste0(exp_gene, " RNA abundance(log2)"))
} else {
  p = p + labs(y=paste0(exp_gene, " ", exp_site, " phosphorylation\n(log2 ratio)"))
  
}
p = p + theme_nogrid()
p = p + theme(title = element_text(size = 18, face = "bold"))
p = p + stat_compare_means(data = tab2p, mapping = aes(x = x, y = y, label = ..p.signif..), symnum.args = symnum.args, ref.group = "Control")     # Add global Anova p-value
p = p + scale_x_discrete(breaks = c("Upstream_Mut_Group", 
                                    "Missense_Group", 
                                    "Trunc_Group",
                                    "Other_Mut_Group",
                                    "Exp_Gene_Mut_Group", 
                                    "Control"),
                         label = c(paste0(mut_gene_upstream, "\nMutated"),
                                   paste0(mut_gene, "\nMissense"),
                                   paste0(mut_gene, "\nTruncation"),
                                   paste0(mut_gene, "\nOther_Mutation"),
                                   paste0(exp_gene, "\nMutated"),
                                   "Control"))
p <- p + guides(fill = F)
if (!is.null(mut_gene_upstream)) {
  p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10, face = "bold"),
                axis.text.x = element_text(size= 15, vjust=0.5, hjust = 0.5, face = "bold"),
                axis.text.y = element_text(colour="black", size=8))
  p
  fn = paste0(subdir3, mut_gene, "_Mut_Detailed_", exp_gene, "_", exp_site, "_with_upstream", paste0(mut_gene_upstream, collapse = "_"), ".pdf")
  pdf(file = fn,  height= 5, width = 3*length(unique(tab2p$group)), useDingbats=FALSE)
  print(p)
  dev.off()
} else {
  p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10, face = "bold"),
                axis.text.x = element_text(size= 15, vjust=0.5, hjust = 0.5, face = "bold"),
                axis.text.y = element_text(colour="black", size=8))
  fn = paste0(subdir3, mut_gene, "_Mut_Detailed_", exp_gene, "_", exp_site, ".pdf")
  pdf(file = fn,  height= 5, width = 3*length(unique(tab2p$group)), useDingbats=FALSE)
  print(p)
  dev.off()
}







# make plot for overall Mut vs WT -----------------------------------------
tab2p <- tab2p_sup
tab2p$x <- tab2p$group
tab2p$x <- factor(tab2p$x, levels = c( "Mut_Group", "Upstream_Mut_Group", "Exp_Gene_Mut_Group", "Control"))
tab2p$y <- as.vector(tab2p$exp_value)
tab2p <- unique(tab2p)
pos <- position_jitter(width = 0.2, seed = 1)
p = ggplot(tab2p, aes(x=x, y=y))
p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.6, size = 2)
p = p + geom_boxplot(aes(fill = group),  color = NA, alpha = 0.4, notch = T)
p = p + scale_fill_manual(values = c("Mut_Group" = set1[1], "Upstream_Mut_Group" = set1[1], "Exp_Gene_Mut_Group" = set1[1], "Control" = "grey50"))
if (exp_site == "PRO") {
  p = p + labs(y=paste0(exp_gene, " protein abundance(log2 ratio)"))
} else if (exp_site == "RNA") {
  p = p + labs(y=paste0(exp_gene, " RNA abundance(log2)"))
} else {
  p = p + labs(y=paste0(exp_gene, " ", exp_site, " phosphorylation\n(log2 ratio)"))
  
}
p = p + theme_nogrid()
p = p + theme(title = element_text(size = 18, face = "bold"))
p = p + stat_compare_means(data = tab2p, mapping = aes(x = x, y = y, label = ..p.signif..), symnum.args = symnum.args, ref.group = "Control")     # Add global Anova p-value
p = p + scale_x_discrete(breaks = c("Upstream_Mut_Group", 
                                    "Mut_Group", 
                                    "Exp_Gene_Mut_Group", 
                                    "Control"),
                         label = c(paste0(mut_gene_upstream, "\nMutated"),
                                   paste0(mut_gene, "\nMutated"),
                                   paste0(exp_gene, "\nMutated"),
                                   "Control"))
p <- p + guides(fill = F)
if (!is.null(mut_gene_upstream)) {
  p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10, face = "bold"),
                axis.text.x = element_text(size= 15, vjust=0.5, hjust = 0.5, face = "bold"),
                axis.text.y = element_text(colour="black", size=8))
  p
  fn = paste0(subdir3, mut_gene, "_", exp_gene, "_", exp_site, "_with_upstream", paste0(mut_gene_upstream, collapse = "_"), ".pdf")
  ggsave(file=fn, height= 3, width = 4, useDingbats=FALSE)
} else {
  p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10, face = "bold"),
                axis.text.x = element_text(size= 20, vjust=0.5, hjust = 0.5, face = "bold"),
                axis.text.y = element_text(colour="black", size=8))
  p
  fn = paste0(subdir3, mut_gene, "_", exp_gene, "_", exp_site, ".pdf")
  pdf(file = fn,  height= 5, width = 3*length(unique(tab2p$group)), useDingbats=FALSE)
  print(p)
  dev.off()
}

tf_tab <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/TF_interactions.txt", data.table = F)
pair_tab <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/analysis_results/dependencies/tables/compile_protein_pair_table/protein_pair_table_v2.txt", data.table = F)
list_downstream <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/list_downstream.txt", data.table = F, col.names = "Gene", header = F)

