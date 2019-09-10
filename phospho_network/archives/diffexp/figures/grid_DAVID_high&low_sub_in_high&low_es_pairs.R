# Yige Wu @ WashU 2018 Jul
# grid view of the same high/low kinase have different high substrates in different pathways

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/pan3can_aes.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# variables ---------------------------------------------------------------
cancers_sort <- c("BRCA", "OV", "CO")
reg_sig <- 0.05
diff_sig <- 0.2
diff_log2fc <- 1
color_cat_man <- c(colors['BRCA'], colors['OV'], colors['COAD'], "#bdbdbd"); names(color_cat_man) <- c("BRCA", "OV", "CO", "NA")

# inputs -------------------------------------------------------------------
sup_cans_tab <- NULL
for (cancer in cancers_sort) {
  sup_tab <- fread(paste0(ppnD, "kinase_activity/tables/fisher_es_pairs/", 
                          cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  sup_tab$Cancer <- cancer
  sup_cans_tab <- rbind(sup_cans_tab, sup_tab)
}
enzyme_type <- "kinase"
sup_cans_tab_pk <- sup_cans_tab[sup_cans_tab$enzyme_type == enzyme_type,]
sup_cans_tab_pk <- data.frame(sup_cans_tab_pk)
sup_cans_tab_pk <- markSigSiteCan(sup_cans_tab_pk, sig_thres = reg_sig, enzyme_type = enzyme_type)
## annotate kinase substrate regulation
sup_cans_tab_pk$regulated <- (sup_cans_tab_pk$coef_sig & sup_cans_tab_pk$fdr_sig)

## annotate up and down regulation
sup_cans_tab_pk$diffexp_type <- "other"
diffexp_sig <- ((!is.na(sup_cans_tab_pk$KSEA_pvalue) & sup_cans_tab_pk$KSEA_pvalue < diff_sig) | (!is.na(sup_cans_tab_pk$diffexp_log2FC & abs(sup_cans_tab_pk$diffexp_log2FC) > diff_log2fc)))
sup_cans_tab_pk$diffexp_type[diffexp_sig & sup_cans_tab_pk$enzyme_direction == "up" & sup_cans_tab_pk$substrate_direction == "up"] <- "up"
sup_cans_tab_pk$diffexp_type[diffexp_sig & sup_cans_tab_pk$enzyme_direction == "down" & sup_cans_tab_pk$substrate_direction == "down"] <- "down"
sum_tab3 <- data.frame(table(unique(sup_cans_tab_pk[, c("pair", "GENE", "diffexp_type", "Cancer")])[, c("GENE","diffexp_type", "Cancer")]))

tab_pk_cans_trans <- sup_cans_tab_pk[sup_cans_tab_pk$SELF == "trans" & !is.na(sup_cans_tab_pk$SELF),]


for (diffexp_type in c("up", "down")) {
  diffexp_type <- "up"
  df <- sum_tab3[sum_tab3$diffexp_type == diffexp_type & sum_tab3$Freq > 1,]
  df_count2 <- data.frame(table(df$GENE))
  df <- df[df$GENE %in% df_count2$Var1[df_count2$Freq == length(cancers_sort)],]
  
  dir1 <- paste0(ppnD, "diffexp/tables/high&low_sub_in_high&low_es_pairs/", diffexp_type, "/")
  paths_cans_genes <- NULL
  # for (enzyme in c("AURKA", "CDK1", "CDK2", "CDK4", "CDK7", "CDK9")) {
    
  for (enzyme in unique(df$GENE)) {
    dir2 <- paste0(dir1, enzyme, "/")

    paths_cans <- NULL
    for (cancer in cancers_sort){
      dir3 <- paste0(dir2, cancer, "/")
      fn <- paste0(dir3, "download.txt")
      if (file.exists(fn)) {
        paths <- fread(input = paste0(dir3, "download.txt"), data.table = F)
        # paths <- paths[paths$Source == "KEGG Pathways",]
        num_genes <- vector(mode = "numeric", length = nrow(paths))
        num_sites <- vector(mode = "numeric", length = nrow(paths))
        desp <- vector(mode = "character", length = nrow(paths))
        
        func_anno_tab <- fread(input = paste0(dir3, "FwDAVID.txt"), data.table = F)
        site_tab <- fread(input = paste0(dir3, diffexp_type, "_", enzyme_type, "_substrate", "_sites_in", cancer, ".txt"), data.table = F, col.names = c("ID", "SUB_MOD_RSD"))
        site_tab <- merge(site_tab, func_anno_tab, all.x = T)
        count <- 0
        for (i in 1:nrow(paths)) {
          count <- count + 1
          p <- paths$Description[i]
          s <- paths$Source[i]
          id <- paths$ID[i]
          
          if (p == "") {
            if (s == "KEGG Pathways") {
              tmp <- as.character(func_anno_tab$KEGG_PATHWAY[grepl(pattern = id, x = as.vector(func_anno_tab$KEGG_PATHWAY))])
              p <- str_split_fixed(string = tmp, pattern = paste0(id, ":"), 2)[,2]
              p <- unique(str_split_fixed(string = p, pattern = ",", 2)[1,1])
            } else if (s == "BioCarta Pathways") {
              tmp <- as.character(func_anno_tab$BIOCARTA[grepl(pattern = id, x = as.vector(func_anno_tab$BIOCARTA))])
              p <- str_split_fixed(string = tmp, pattern = paste0(id, ":"), 2)[,2]
              p <- unique(str_split_fixed(string = p, pattern = ",", 2)[1,1])
            }
          }
          if (s == "KEGG Pathways") {
            num_genes[count] <- nrow(func_anno_tab[grepl(pattern = p, x = as.vector(func_anno_tab$KEGG_PATHWAY)),])
            num_sites[count] <- nrow(site_tab[grepl(pattern = p, x = as.vector(site_tab$KEGG_PATHWAY)),])
          } else if (s == "BioCarta Pathways") {
            num_genes[count] <- nrow(func_anno_tab[grepl(pattern = p, x = as.vector(func_anno_tab$BIOCARTA)),])
            num_sites[count] <- nrow(site_tab[grepl(pattern = p, x = as.vector(site_tab$BIOCARTA)),])
          }
          desp[count] <- p
        }
        paths$num_genes <- num_genes
        paths$num_sites <- num_sites
        paths$Description <- desp
        paths$Cancer <- cancer
        paths_cans <- rbind(paths_cans, paths)
        
      }
    }
    paths_cans$GENE <- enzyme
    paths_cans_genes <- rbind(paths_cans_genes, paths_cans)
  }
  paths_cans_genes <- paths_cans_genes[paths_cans_genes$num_sites >= 4,]
  for (s in unique(paths_cans_genes$Source)) {
    tab2p <- data.frame(table(paths_cans_genes[paths_cans_genes$Source == s, c("Cancer", "Description", "GENE")]))
    tab2p$x_print <- paste0(tab2p$Description, "_", tab2p$Cancer)
    tab2p$fill <- tab2p$Cancer
    tab2p$fill[tab2p$Freq == 0] <- NA
    
    p <- ggplot()
    p <- p + geom_tile(data = tab2p, mapping = aes(x = x_print, y = GENE, fill = fill))
    p <- p + scale_fill_manual(values = color_cat_man)
    p <- p + facet_grid(.~Description, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p <- p + theme_nogrid()
    p <- p + theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
                   strip.text.x = element_text(size = 5, angle = 90),
                   panel.spacing.x = unit(0, "lines"))
    p 
    fn = paste0(makeOutDir(resultD = resultD),'3can_', diffexp_type,'_kinase_substrates_', s, '_across_cancers.pdf')
    ggsave(file=fn, height=5, width=10)
  }
}
