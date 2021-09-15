# Yige Wu @ WashU 2020 Nov.

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/"
setwd(dir_base)
source("./phospho-signaling_analysis/load_pkgs.R")
source("./phospho-signaling_analysis/functions.R")
source("./phospho-signaling_analysis/plotting_dependencies.R")
source("./phospho-signaling_analysis/dependencies/tables/load_TCGA_pathways.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# set variables -----------------------------------------------------------
sig_fdr_thres <- 0.05

# input dependencies ------------------------------------------------------
result_df <- fread(input = "./analysis_results/phospho_network/genoalt/tables/test_PDA_SMG_mut_impact_global/20201130.v1/PDAC_not_silent_mut_impact_PRO.tsv", data.table = F)

# make plot data frame ----------------------------------------------------
result_df <- result_df %>%
  mutate(is_significant = (fdr <= sig_fdr_thres))
  
## filter for affected genes to show
genes_affected_plot <- result_df$SUB_GENE[result_df$is_significant]
plot_data_df <- result_df %>%
  filter(SUB_GENE %in% genes_affected_plot) %>%
  mutate(log10_FDR = (-log10(fdr)))
## annotate affected genes to pathways
genes2pathways <- map2TCGApathwaways(gene_list = unique(plot_data_df$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
genes2pathways_df <- data.frame(SUB_GENE = rep(x = names(genes2pathways), sapply(X = genes2pathways, FUN = function(x) length(x))), 
                                SUB_GENE.path = unlist(genes2pathways, use.names = F))
### manually remove some annotation
genes2pathways_filtered_df <- genes2pathways_df %>%
  arrange(SUB_GENE) %>%
  filter(SUB_GENE == "BCL2L1" & SUB_GENE.path == "RTK RAS" | SUB_GENE != "BCL2L1") %>%
  filter(SUB_GENE == "EPHA2" & SUB_GENE.path == "RTK RAS" | SUB_GENE != "EPHA2") %>%
  filter(SUB_GENE == "MET" & SUB_GENE.path == "RTK RAS" | SUB_GENE != "MET") %>%
  filter(SUB_GENE == "PCNA" & SUB_GENE.path == "Cell Cycle" | SUB_GENE != "PCNA") %>%
  filter(SUB_GENE == "SFN" & SUB_GENE.path == "TP53" | SUB_GENE != "SFN") %>%
  filter(SUB_GENE == "RBX1" & SUB_GENE.path == "Cell Cycle" | SUB_GENE != "RBX1")
## map pathway
subgene_path_vec <- mapvalues(x = plot_data_df$SUB_GENE, from = genes2pathways_filtered_df$SUB_GENE, to = as.vector(genes2pathways_filtered_df$SUB_GENE.path))
subgene_path_vec[subgene_path_vec %in% plot_data_df$SUB_GENE] <- "Other"
plot_data_df$SUB_GENE.path <- factor(x = subgene_path_vec, levels = c("RTK RAS", "Cell Cycle", "Mismatch repair", "PI3K", "TP53", "HIPPO", "WNT", "Other"))

## adjust median difference
summary(plot_data_df$meddiff)
cap <- min(max(abs(plot_data_df$meddiff), na.rm = T), 1)
plot_data_df <- plot_data_df %>%
  mutate(meddiff_capped = ifelse(meddiff > cap | meddiff < (-meddiff), cap, meddiff))

# make colors -------------------------------------------------------------
sig_colors <- c("black", "white")
names(sig_colors) <- c(paste0('FDR<', sig_fdr_thres), paste0('FDR>', sig_fdr_thres))
colors_exp <- RdBu1024

# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plot_data_df, mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_FDR, shape = SELF),
                                                       alpha = 0.6, color = "black")
p <- p  + scale_shape_manual(values = c("cis" = 22, "trans" = 21))
p <- p + scale_color_manual(values = sig_colors)
p = p + scale_fill_gradientn(name= "protein\nchange", na.value=NA, colours=colors_exp, limit=c(-cap,cap))
p <- p + facet_grid(. ~ SUB_GENE.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
p <- p + xlab("Affected Protein") + ylab("Mutated Gene")
p <- p + guides(colour = guide_legend(title = "FDR"))
p <- p + theme_bw()
p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
p <- p + theme(strip.text.x = element_text(angle = 90))
p <- p + theme(legend.position = "bottom")
p <- p + guides(size = guide_legend(nrow = 3, byrow = T))
p

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "PDA_Mut_Impact_PRO.png")
png(file2write, width = 1000, height = 600, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "PDA_Mut_Impact_PRO.pdf")
pdf(file2write, width = 6, height = 4, useDingbats = F)
print(p)
dev.off()
