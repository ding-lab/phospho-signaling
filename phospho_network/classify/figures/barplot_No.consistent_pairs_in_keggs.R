# Yige Wu @ WashU 2018 Apr
# showing # of consistent pairs in hallmark pathways


# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/pan3can_aes.R')


# inputs ------------------------------------------------------------------
load("~/Box Sync/pan3can_shared_data/analysis_results/2015-08-01_Gene_Set.RData")
keywords <- c("cell cycle", "Hippo", "Notch", "PI3K", "TGF", "Ras", "p53", "Wnt")


# get kegg pathways -------------------------------------------------------
all_keggs <- names(KEGG)
keggs <- NULL
for (key in keywords) {
  keggs <- c(keggs, all_keggs[grepl(key, all_keggs, ignore.case = T)])
}
keggs <- unique(keggs)

#for (cancer in c("BRCA")) {
for (cancer in c("BRCA", "OV", "CO")) {
  ksea_diffexp <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_regression/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  ## get consistently high/low pairs
  cons_tab <- ksea_diffexp[!is.na(ksea_diffexp$substrate_direction) & !is.na(ksea_diffexp$consistent) & ksea_diffexp$consistent,]
  cons_up_tab <- cons_tab[cons_tab$substrate_direction == "up",]
  cons_down_tab <- cons_tab[cons_tab$substrate_direction == "down",]
  
  ##get the genes in those pairs
  cons_all_genes <- unique(c(as.vector(cons_tab$GENE), as.vector(cons_tab$SUB_GENE)))
  cons_up_genes <- unique(c(as.vector(cons_up_tab$GENE), as.vector(cons_up_tab$SUB_GENE)))
  cons_down_genes <- unique(c(as.vector(cons_down_tab$GENE), as.vector(cons_down_tab$SUB_GENE)))
  
  ## initiate the table to plot
  kegg_num <- data.frame(kegg_name = keggs,
                         kegg_num_genes = sapply(keggs, FUN = function(k, l) length(l[[k]]), l = KEGG))
  kegg_num_genes_cons <- vector(mode = "numeric", length = length(keggs)); names(kegg_num_genes_cons) <- keggs
  kegg_num_genes_cons_up <- vector(mode = "numeric", length = length(keggs)); names(kegg_num_genes_cons_up) <- keggs
  kegg_num_genes_cons_down <- vector(mode = "numeric", length = length(keggs)); names(kegg_num_genes_cons_down) <- keggs
  
  for (kegg_tmp in keggs) {
    cons_all_genes_tmp <- cons_all_genes[cons_all_genes %in% KEGG[[kegg_tmp]]]
    kegg_num_genes_cons[kegg_tmp] <- length(cons_all_genes_tmp)
    cons_up_genes_tmp <- cons_up_genes[cons_up_genes %in% KEGG[[kegg_tmp]]]
    kegg_num_genes_cons_up[kegg_tmp] <- length(cons_up_genes_tmp)
    cons_down_genes_tmp <- cons_down_genes[cons_down_genes %in% KEGG[[kegg_tmp]]]
    kegg_num_genes_cons_down[kegg_tmp] <- length(cons_down_genes_tmp)
    if (length(cons_all_genes_tmp) >0) {
      print(paste0(cancer, ":", kegg_tmp, ":", length(cons_all_genes_tmp)))
      cons_tab[, kegg_tmp] <- ifelse(cons_tab$GENE %in% cons_all_genes_tmp | cons_tab$SUB_GENE %in% cons_all_genes_tmp, TRUE, FALSE)
    }
  }
  kegg_num <- cbind(kegg_num, data.frame(kegg_num_genes_cons = kegg_num_genes_cons, kegg_num_genes_cons_up = kegg_num_genes_cons_up, kegg_num_genes_cons_down = kegg_num_genes_cons_down))
  kegg_num$up_gene_ratio <- kegg_num$kegg_num_genes_cons_up/kegg_num$kegg_num_genes_cons_down
  cons_tab_m <- rbind(data.frame(kegg_name = keggs,
                                 substrate_direction = "up",
                                 num_genes = kegg_num_genes_cons_up,
                                 up_gene_ratio = kegg_num$up_gene_ratio),
                      data.frame(kegg_name = keggs,
                                 substrate_direction = "down",
                                 num_genes = -(kegg_num_genes_cons_down),
                                 up_gene_ratio = kegg_num$up_gene_ratio))
  cons_tab_m$kegg_name <- factor(cons_tab_m$kegg_name, levels = as.vector(kegg_num$kegg_name)[order(kegg_num$up_gene_ratio)])

  p <- ggplot()
  p <- p + geom_bar(data=cons_tab_m, aes(y = num_genes, x = kegg_name, fill = substrate_direction), stat="identity")
  p <- p + geom_text(data=cons_tab_m, aes(y = ifelse(up_gene_ratio > 1, 10, -10), x = kegg_name, label = kegg_name), size = 4)
  p <- p + scale_fill_manual(values = c("down" = "#1F78B4", "up" = "#E31A1C"))
  p <- p + theme_bw() + theme_nogrid()
  p <- p + coord_flip()
  p <- p + xlab("enzyme")+ylab("number of genes in consistently up/down enzyme-substrate pairs")
  p <- p + theme(axis.title=element_text(size=10))
  p <- p + theme(axis.text.x = element_text(colour="black", size = 10),
                 axis.text.y = element_blank())
  p
  ggsave(filename = paste0(makeOutDir(), "No.consistent_highorlow_pairs_within_keggs_in_", cancer, ".pdf"), 
         height = 6, width = 6)
}

