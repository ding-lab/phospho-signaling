# Yige Wu @ WashU 2021 Mar

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/"
setwd(dir_base)
source("./phospho-signaling_analysis/load_pkgs.R")
source("./phospho-signaling_analysis/functions.R")
source("./phospho-signaling_analysis/plotting_dependencies.R")
library(ggpubr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the mutation data
maf_df <- fread(data.table = F, input = "~/Box/Data_Freeze_1.0/CPTAC_ccRCC_Combined/Somatic_mutation/CPTAC_ccRCC_combined_somatic_mutation.maf_v1.0.tsv")
## input protein data
pro_df <- fread(data.table = F, input = "~/Box/Data_Freeze_1.0/CPTAC_ccRCC_Combined/Proteome/UMich_pipeline/ccRCC_DIA_Combined_norm_03042021.tsv", fill = T)
## input test results
result_df <- fread(data.table = F, input = "./analysis_results/phospho_network/genoalt/tables/test_ccRCC_SMG_mut_impact_global_protein_RNA/20210331.v1/ccRCC_not_silent_mut_impact_RNA.tsv")
# result_df <- fread(data.table = F, input = "./analysis_results/phospho_network/genoalt/tables/test_ccRCC_SMG_mut_impact_global_protein_RNA/20210331.v1/ccRCC_truncation_mut_impact_PRO.tsv")
# result_df <- fread(data.table = F, input = "./analysis_results/phospho_network/genoalt/tables/test_ccRCC_SMG_mut_impact_global_protein_RNA/20210331.v1/ccRCC_not_silent_mut_impact_PRO.tsv")

# set parameters for plottting --------------------------------------------------
## make plotting parameters
pos <- position_jitter(width = 0.2, seed = 1)
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
colors_samplegroup <- RColorBrewer::brewer.pal(n = 5, name = "Set1")
names(colors_samplegroup) <- c("Missense\nTumor", "Truncation\nTumor", "NAT", "WT\nTumor", "Other Mut\nTumor")
## filter
toplot_df <- result_df %>%
  # filter(fdr < 0.05) %>%
  # filter(fdr_by_gene < 0.05) %>%
  filter(SELF == "cis") %>%
  select(GENE, SUB_GENE)
## log 2 protein values
exp_head <- data.frame(gene_name = pro_df$Genes)
colnames_pro <- colnames(pro_df)
exp_data <- pro_df[, colnames_pro[grepl(pattern = "_T|_N", x = colnames_pro) & !grepl(pattern = "\\.1", x = colnames_pro)]]
exp_data <- log2(exp_data)
exp_df <- cbind(exp_head, exp_data)

# plot by row above -------------------------------------------------------
i <- 1
for (i in 1:nrow(toplot_df)) {
  gene_mut <- toplot_df$GENE[i]
  gene_exp <- toplot_df$SUB_GENE[i]
  
  ## get mutation status
  case2mut_wide_df <- get_mutation_class_sim_matrix(pair_tab = gene_mut, maf = maf_df)
  case2mut_long_df <- melt(data = case2mut_wide_df, id.var = c("Hugo_Symbol"))
  case2mut_long_df <- case2mut_long_df %>%
    mutate(variant_class = ifelse(value %in% c("Missense", "Truncatin"), value, "Other Mut")) %>%
    mutate(variant_class = ifelse(grepl(x = value, pattern = "Truncation", ignore.case = T), "Truncation", variant_class)) %>%
    mutate(sample_id = paste0(variable, "_T"))
  ## get expression data
  plot_data_wide_df <- exp_df %>%
    dplyr::filter(gene_name == gene_exp)
  plot_data_long_df <- melt(data = plot_data_wide_df)
  ## merge
  plot_data_df <- merge(x = plot_data_long_df, 
                        y = case2mut_long_df, by.x = c("variable"), by.y = c("sample_id"), all.x = T)
  plot_data_df <- plot_data_df %>%
    mutate(sample_group = ifelse(!is.na(variant_class), paste0(variant_class, "\nTumor"),
                                  ifelse(grepl(x = variable, pattern = "_T"), "WT\nTumor", "NAT"))) %>%
    mutate(y_plot = value.x) %>%
    mutate(x_plot = sample_group)
  plot_data_df$x_plot <- factor(x = plot_data_df$x_plot, levels = c("NAT", "WT\nTumor", "Truncation\nTumor", "Missense\nTumor", "Other Mut\nTumor"))

  ## plot
  p <- ggplot(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 0.6)
  p = p + geom_boxplot(width=.1, alpha = 0.8)
  p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.7, size = 1.5)
  p = p + stat_compare_means(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot, label = ..p.signif..), 
                             symnum.args = symnum.args, ref.group = "NAT")     # Add global Anova p-value
  p <- p + scale_fill_manual(values = colors_samplegroup)
  p <- p + ggtitle(paste0(gene_mut, " mutation status", " ~ " , gene_exp, " protein level"))
  p <- p + ylab(paste0(gene_exp, " protein level (log2)"))
  p <- p + theme_classic(base_size = 15)
  p <- p + theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size = 12))
  # p
  file2write <- paste0(dir_out, gene_mut, "_Mut.", gene_exp, "_Protein", ".png")
  png(file2write, width = 700, height = 600, res = 150)
  print(p)
  dev.off()
}
