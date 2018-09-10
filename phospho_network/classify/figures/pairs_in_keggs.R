

source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/pan3can_aes.R')

upstream_count <- data.frame(table(us_diffexp_nona[us_diffexp_nona$SELF == "trans", c("SUBSTRATE", "upstream")]))
upstream_count <- upstream_count[upstream_count$Freq>0,]
upstream_type_count <- data.frame(table(upstream_count$SUBSTRATE))
upstream_type_duo <- upstream_type_count[upstream_type_count$Freq>=2,]

for (substrate in upstream_type_duo$Var1) {
  tmp <- us_diffexp_nona[us_diffexp_nona$SUBSTRATE == substrate & any(us_diffexp_nona$consistent),]
  print(tmp)
}

us_diffexp_nona_trans <- us_diffexp_nona[us_diffexp_nona$SELF == "trans",]
for (cancer in c("BRCA", "OV", "CO")) {
  us_all_genes <- unique(c(as.vector(us_diffexp_nona_trans$KINASE[us_diffexp_nona_trans$Cancer == cancer]), as.vector(us_diffexp_nona_trans$SUBSTRATE[us_diffexp_nona_trans$Cancer == cancer])))
  kegg_anno <- data.frame(gene = us_all_genes)
  for (kegg_tmp in keggs) {
    gene_tmp <- us_all_genes[us_all_genes %in% KEGG[[kegg_tmp]]]
    if (length(gene_tmp) >0) {
      print(paste0(cancer, ":", kegg_tmp, ":", length(gene_tmp)))
      kegg_anno[, kegg_tmp] <- ifelse(kegg_anno$gene %in% gene_tmp, TRUE, FALSE)
    }
  }
  kegg_anno <- kegg_anno[rowSums(kegg_anno == TRUE) > 0,]
  kegg_anno_m <- melt(kegg_anno)
  kegg_anno_m <- kegg_anno_m[kegg_anno_m$value == TRUE,]
  us_tmp <- unique(us_diffexp_nona_trans[us_diffexp_nona_trans$Cancer == cancer , c("KINASE", "SUBSTRATE", "upstream_direction", "downstream_direction", "consistent")])
  us_tmp <- merge(us_tmp, kegg_anno_m, by.x = c("KINASE"), by.y = c("gene"), all = T)
  colnames(us_tmp)[colnames(us_tmp) %in% c("variable", "value")] <- c("kegg_kin", "iskegg_kin")
  us_tmp <- merge(us_tmp, kegg_anno_m, by.x = c("SUBSTRATE"), by.y = c("gene"), all = T)
  colnames(us_tmp)[colnames(us_tmp) %in% c("variable", "value")] <- c("kegg_sub", "iskegg_sub")
  us_tmp <- us_tmp[(!is.na(us_tmp$iskegg_kin) & !is.na(us_tmp$iskegg_sub)),]
  us_tmp <- us_tmp[us_tmp$iskegg_kin == TRUE & us_tmp$iskegg_sub == TRUE,]
  p <- ggplot(us_tmp)
  p <- p + geom_tile(aes(x = KINASE, y = SUBSTRATE, fill = downstream_direction, color = upstream_direction, alpha = 0.5*((consistent == TRUE) + 1)))
  p <- p + facet_grid(kegg_kin~kegg_sub, drop = T, space = "free",scales = "free")
  p <- p + theme_bw()
  p = p + theme(strip.text.y = element_text(size = 7, angle = 0), strip.text.x = element_text(size = 5))
  p = p + theme(axis.text.x = element_text(size = 5, angle = -30, hjust = 0, vjust = 0.1))
  p <- p + scale_fill_manual(values = c("up" = set1[1], "down" = set1[2]))
  p <- p + scale_color_manual(values = c("up" = set1[1], "down" = set1[2]))
  p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
  p
  ggsave(filename = paste0(resultDnow, "ks_pairs_across_keggs_in_", cancer, ".pdf"), 
         height = 12, width = 12)
  
  us_tmp2 <- unique(us_diffexp_nona_trans[us_diffexp_nona_trans$Cancer == cancer , c("KINASE", "SUBSTRATE", "SUB_MOD_RSD", "upstream_direction", "downstream_direction", "consistent")])
  us_tmp2 <- merge(us_tmp2, kegg_anno_m, by.x = c("KINASE"), by.y = c("gene"), all = T)
  colnames(us_tmp2)[colnames(us_tmp2) %in% c("variable", "value")] <- c("kegg_kin", "iskegg_kin")
  us_tmp2 <- merge(us_tmp2, kegg_anno_m, by.x = c("SUBSTRATE"), by.y = c("gene"), all = T)
  colnames(us_tmp2)[colnames(us_tmp2) %in% c("variable", "value")] <- c("kegg_sub", "iskegg_sub")
  us_tmp2 <- us_tmp2[(!is.na(us_tmp2$iskegg_kin) & !is.na(us_tmp2$iskegg_sub)),]
  us_tmp2 <- us_tmp2[us_tmp2$iskegg_kin == TRUE & us_tmp2$iskegg_sub == TRUE,]
  us_tmp2 <- us_tmp2[us_tmp2$kegg_kin == as.vector(us_tmp2$kegg_sub),]
  us_tmp2$phosphosite <- paste0(us_tmp2$SUBSTRATE, ":", us_tmp2$SUB_MOD_RSD)
  us_tmp2$id <- paste0(us_tmp2$KINASE, us_tmp2$SUBSTRATE, us_tmp2$consistent)
  keeprow <- (us_tmp2$consistent | (!us_tmp2$consistent & !duplicated(us_tmp2$id)) )
  
  p <- ggplot(us_tmp2[keeprow,])
  p <- p + geom_tile(aes(x = KINASE, y = phosphosite, fill = downstream_direction, color = upstream_direction, alpha = 0.5*((consistent == TRUE) + 1)))
  p <- p + facet_grid(kegg_kin~kegg_sub, drop = T, space = "free",scales = "free")
  p <- p + theme_bw()
  p = p + theme(strip.text.y = element_text(size = 7, angle = 0), strip.text.x = element_text(size = 5))
  p = p + theme(axis.text.x = element_text(size = 8, angle = -30, hjust = 0, vjust = 0))
  p <- p + scale_fill_manual(values = c("up" = set1[1], "down" = set1[2]))
  p <- p + scale_color_manual(values = c("up" = set1[1], "down" = set1[2]))
  p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
  p
  ggsave(filename = paste0(resultDnow, "ks_pairs_within_keggs_in_", cancer, ".pdf"), 
         height = 12, width = 12)
}

