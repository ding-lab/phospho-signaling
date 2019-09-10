# Yige Wu @ WashU 2018 Jul
# 

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(igraph)

# variables ---------------------------------------------------------------
cancers_sort <- c("BRCA", "OV", "CO")
reg_sig <- 0.05
diff_sig <- 0.2
diff_log2fc <- 1

# inputs -------------------------------------------------------------------
sup_cans_tab <- NULL
for (cancer in cancers_sort) {
  sup_tab <- fread(paste0(ppnD, "kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/", 
                          cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  sup_tab$Cancer <- cancer
  sup_cans_tab <- rbind(sup_cans_tab, sup_tab)
}

# classify ----------------------------------------------------------------
sup_cans_tab_pk <- sup_cans_tab[sup_cans_tab$enzyme_type == "kinase",]
sup_cans_tab_pk <- markSigSiteCan(sup_cans_tab_pk, sig_thres = reg_sig, enzyme_type = "kinase")
sup_cans_tab_pp <- sup_cans_tab[sup_cans_tab$enzyme_type == "phosphotase",]
sup_cans_tab_pp <- markSigSiteCan(sup_cans_tab_pp, sig_thres = reg_sig, enzyme_type = "phosphotase")
sup_cans_tab <- rbind(sup_cans_tab_pk, sup_cans_tab_pp)
## annotate kinase substrate regulation
sup_cans_tab$regulated <- (sup_cans_tab$coef_sig & sup_cans_tab$fdr_sig)
## annotate up and down regulation
sup_cans_tab$diffexp_type <- "other"
diffexp_sig <- ((!is.na(sup_cans_tab$KSEA_pvalue) & sup_cans_tab$KSEA_pvalue < diff_sig) | (!is.na(sup_cans_tab$diffexp_log2FC & abs(sup_cans_tab$diffexp_log2FC) > diff_log2fc)))
sup_cans_tab$diffexp_type[diffexp_sig & sup_cans_tab$enzyme_direction == "up" & sup_cans_tab$substrate_direction == "up"] <- "up"
sup_cans_tab$diffexp_type[diffexp_sig & sup_cans_tab$enzyme_direction == "down" & sup_cans_tab$substrate_direction == "down"] <- "down"
sup_cans_tab_pk <- sup_cans_tab[sup_cans_tab$enzyme_type == "kinase",]
sup_cans_tab_pp <- sup_cans_tab[sup_cans_tab$enzyme_type == "phosphotase",]

graph_stat_cans <- NULL
sup_cans_tab_pk_trans <- sup_cans_tab_pk[sup_cans_tab_pk$SELF == "trans" & !is.na(sup_cans_tab_pk$SELF),]
for (cancer in unique(sup_cans_tab_pk_trans$Cancer)) {
  print(cancer)
  tab4igraph <- sup_cans_tab_pk_trans
  tab4igraph <- tab4igraph[tab4igraph$regulated & tab4igraph$Cancer == cancer,]
  tab4igraph <- unique(tab4igraph[,c("pair", "GENE", "SUB_GENE", "SUB_MOD_RSD")])
  relations <- data.frame(from = as.vector(tab4igraph$GENE), to = as.vector(tab4igraph$SUB_GENE))
  relations <- relations[!is.na(relations$from) & !is.na(relations$to),]
  nodes <- data.frame(name = unique(c(as.vector(relations$from), as.vector(relations$to))))
  nodes <- nodes[!is.na(nodes$name),]
  nodes <- data.frame(name = nodes)
  
  g <- graph_from_data_frame(relations, directed = TRUE, vertices = nodes)
  # print(g, e=TRUE, v = T)
  
  ndbtw <- betweenness(g, normalized = T)
  ndbtw <- ndbtw[order(ndbtw, decreasing = T)]
  print(head(ndbtw))
  tmp <- data.frame(name = names(ndbtw), betweenness = ndbtw)
  graph_stat <- tmp
  
  ndcls_in <- closeness(g, mode="in", normalized = T)
  ndcls_in <- ndcls_in[order(ndcls_in, decreasing = T)]
  tmp <- data.frame(name = names(ndcls_in), closeness_in = ndcls_in)
  graph_stat <- merge(graph_stat, tmp, all = T)
  
  ndcls_out <- closeness(g, mode="out", normalized = T)
  ndcls_out <- ndcls_out[order(ndcls_out, decreasing = T)]
  tmp <- data.frame(name = names(ndcls_out), closeness_out = ndcls_out)
  graph_stat <- merge(graph_stat, tmp, all = T)
  
  ndcls_total <- closeness(g, mode="total", normalized = T)
  ndcls_total <- ndcls_out[order(ndcls_total, decreasing = T)]
  tmp <- data.frame(name = names(ndcls_total), closeness_total = ndcls_total)
  graph_stat <- merge(graph_stat, tmp, all = T)
  
  nddg_in <- degree(g, mode = "in")
  tmp <- data.frame(name = names(nddg_in), degree_in = nddg_in)
  graph_stat <- merge(graph_stat, tmp, all = T)
  
  nddg_out <- degree(g, mode = "out")
  tmp <- data.frame(name = names(nddg_out), degree_out = nddg_out)
  graph_stat <- merge(graph_stat, tmp, all = T)
  
  nddg_total <- degree(g, mode = "total")
  tmp <- data.frame(name = names(nddg_total), degree_total = nddg_total)
  graph_stat <- merge(graph_stat, tmp, all = T)

  graph_stat$Cancer <- cancer
  graph_stat_cans <- rbind(graph_stat_cans, graph_stat)
}
write.table(x = graph_stat_cans, file = paste0(makeOutDir(resultD = resultD), "igraph_stats.txt"), quote = F, row.names = F)
