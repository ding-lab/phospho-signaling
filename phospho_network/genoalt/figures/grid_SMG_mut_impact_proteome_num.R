# Yige Wu @ WashU 2018 Feb
# draw a grid showing the overview of the mutational impact of SMGs

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
library(pheatmap)


# input pair annotation table ---------------------------------------------
pair_tab_annotated <- fread(input = "./cptac2p/analysis_results/dependencies/tables/compile_protein_pair_table/protein_pair_table.txt", data.table = F)
colnames(pair_tab_annotated)

# set variable ------------------------------------------------------------
pair_cat <- colnames(pair_tab_annotated)[!(colnames(pair_tab_annotated) %in% c("GENE", "SUB_GENE", "pair_pro"))]
pair_cat
pair_cat_signaling <- c("SUB_GENE.is_TF_downstream", "SUB_GENE.is_kinase_substrate", "SUB_GENE.is_phosphatase_substrate", "SUB_GENE.is_complex_partner")
## choose the type of mutation-affected protein to draw
## the columns to plot as row annotation for the heatmap
# cat2p <- pair_cat
cat2p <- pair_cat_signaling
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")

# input mutation impact table ---------------------------------------------
mut_impact_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/genoalt/tables/test_SMG_mut_impact_proteome/SMG_mut_impact_tab.txt", data.table = F)
mut_impact_tab$cat2p_sum <- rowSums(mut_impact_tab[, cat2p])
mut_impact_tab$GENE.is_SMG <- get_SMG_by_cancer(gene_vector = mut_impact_tab$GENE, cancer_vector = mut_impact_tab$cancer)


# construct heatmap value matrix ------------------------------------------
## the number of proteins affected by SMGs
mut_impact_tab2p <- mut_impact_tab %>%
  filter(fdr_by_gene < 0.05) %>%
  filter(cat2p_sum > 0 | SELF == "cis") %>%
  filter(GENE.is_SMG == T) %>%
  filter(num >= 5) %>%
  filter(cancer %in% cancers2process)
mut_impact_tab2p$SUB_GENE_cat <- sapply(X = 1:nrow(mut_impact_tab2p), FUN = function(i, df, col_names) paste0(unlist(df[i, col_names]), collapse = "_"), df = mut_impact_tab2p, col_names = cat2p)
mut_impact_tab2p$SUB_GENE_cat[mut_impact_tab2p$SELF == "cis"] <- sort(unique(mut_impact_tab2p$SUB_GENE_cat))[1]
mut_impact_tab2p[mut_impact_tab2p$SELF == "cis", cat2p] <- FALSE
df1 <- data.frame(table(unique(mut_impact_tab2p[, c("affected_exp_type", "GENE", "cancer", cat2p, "SUB_GENE")])[, c("affected_exp_type", "GENE", "cancer", cat2p)]))
df1$SUB_GENE_cat <- sapply(X = 1:nrow(df1), FUN = function(i, df, col_names) paste0(unlist(df[i, col_names]), collapse = "_"), df = df1, col_names = cat2p)

df1$GENE.is_SMG <- get_SMG_by_cancer(gene_vector = df1$GENE, cancer_vector = df1$cancer)
df1 <- df1 %>%
  mutate(GENE_cancer = paste0(GENE, "_", cancer)) %>%
  filter(GENE.is_SMG == T)

df1$row_id <- paste0(df1$SUB_GENE_cat, "_", df1$affected_exp_type)
df1 %>%
  head()
df2 <- dcast(data = df1, row_id ~  GENE_cancer, value.var = "Freq")
rownames(df2) <- df2$row_id

# construct heat map row annotation ---------------------------------------
## a matrix of the relation of the affected protein to the mutated protein
## add the expression type
row_anno <- unique(df1[, c("row_id", "affected_exp_type", cat2p)])
# row_anno <- merge(row_anno, unique(mut_impact_tab2p[, c("SUB_GENE_cat", cat2p)]), by = c("SUB_GENE_cat"), all.x = T)
row_anno$SUB_GENE_cat <- NULL
row_anno %>%
  head()
row_anno <- row_anno[order(row_anno$affected_exp_type, decreasing = T),]
rownames(row_anno) <- row_anno$row_id
row_order <- row_anno$row_id
row_anno$row_id <- NULL
row_anno$SUB_GENE_cat <- NULL

# construct heatmap column annotation -------------------------------------
## add the SMG cancer type
col_anno <- unique(df1[, c("GENE_cancer", "cancer")])
col_anno %>%
  head()
rownames(col_anno) <- col_anno$GENE_cancer
col_anno <- col_anno[order(col_anno$cancer),]
col_order <- col_anno$GENE_cancer
col_anno$GENE_cancer <- NULL

# sorting and filtering matrix ----------------------------------------------------------
col_order
row_order
mat_value <- df2[row_order, col_order]

row_keep <- rownames(mat_value)[rowSums(mat_value) > 0]
col_keep <- colnames(mat_value)[colSums(mat_value) > 0]

mat_value_keep <- mat_value[row_keep, col_keep]
col_anno_keep <- col_anno[col_keep, "cancer"] 
col_anno_keep <- data.frame(cancer = col_anno_keep); rownames(col_anno_keep) <- col_keep
row_anno_keep <- row_anno[row_keep,]

mat_label <- mat_value_keep
mat_label[mat_label == 0] <- ""
head(mat_label)


# annotation colors -------------------------------------------------------
row_anno_colors <- list()
for (row_col in cat2p) {
  row_anno_colors[[row_col]] <- c("TRUE" = "black", "FALSE" = "white")
}
tmp <- rainbow(length(unique(col_anno_keep$cancer)))
names(tmp) <- unique(col_anno_keep$cancer)
row_anno_colors[["cancer"]] <- tmp

# draw heatmap ------------------------------------------------------------
my_heatmap <- pheatmap(mat_value_keep, 
                       color = c("white", col_paletteR(100)),
                       na_col = "white", 
                       annotation_col = col_anno_keep,
                       annotation_row = row_anno_keep, 
                       annotation_colors = row_anno_colors,
                       display_numbers = mat_label,
                       cellwidth = 10,
                       cellheight = 8,
                       cluster_rows=F, cluster_cols = F, border_color = "grey",
                       show_rownames = F,
                       show_colnames = T)
fn = paste(makeOutDir(resultD = resultD), "SMG_mut_impact_proteome_numbers.pdf",sep = "")

save_pheatmap_pdf(x = my_heatmap, 
                  filename = fn, 
                  width = 12, height = 7)

