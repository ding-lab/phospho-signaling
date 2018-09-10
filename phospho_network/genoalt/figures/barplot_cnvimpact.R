# Yige Wu @ WashU 2018 July
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
sig_thres <- 0.05
top_num <- 15
networkin_thres <- 5

# inputs ------------------------------------------------------------------
## input enzyme-substrate table
# ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD_extended.csv"))
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor/omnipath_networkin_DEPOD_SignorNotSiteMapped.csv"))
# ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source == "NetKIN" | (ptms_site_pairs_sup$Source != "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= networkin_thres),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)


# count per enzyme  ------------------------
count_enzyme_cans <- NULL
for (cancer in cancers_sort) {
  mut_cnv_tab <- fread(input = paste0(ppnD, "genoalt/tables/mut_cnv_impact/", cancer, "_mut_cnv_tab.txt"), data.table = F)
  mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE)
  mut_cnv_tab <- mut_cnv_tab[mut_cnv_tab$pair_pro %in% ptms_site_pairs_sup$pair_pro,]
  mut_cnv_tab$pair <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE, ":", mut_cnv_tab$SUB_MOD_RSD)
  mut_cnv_tab <- merge(mut_cnv_tab, unique(ptms_site_pairs_sup[, c("GENE", "enzyme_type")]), all.x = T)
  
  ##filter out bad annotated enzyme type
  mut_cnv_tab <- mut_cnv_tab[!(mut_cnv_tab$GENE == "PTPRG" & mut_cnv_tab$enzyme_type == "kinase"),]
  
  ## filter mutation-impacted kinase-substrate pairs
  mut_cnv_mut_tab <- mut_cnv_tab
  mut_cnv_mut_tab <- mut_cnv_tab[mut_cnv_tab$p_mut > 0 & mut_cnv_tab$p_mut < sig_thres,]
  mut_cnv_mut_tab <- mut_cnv_mut_tab[(mut_cnv_mut_tab$meddiff_bottom_mut > 0) == (mut_cnv_mut_tab$meddiff_upper_mut > 0),]
  mut_cnv_mut_tab <- mut_cnv_mut_tab[order(mut_cnv_mut_tab$p_mut),]
  count_enzyme_mut <- data.frame(table(mut_cnv_mut_tab$GENE))
  if (nrow(count_enzyme_mut) > 0){
    count_enzyme_mut$genoalt_type <- "mutation"
  }
  
  ## filter amplification-impacted kinase-substrate pairs
  mut_cnv_amp_tab <- mut_cnv_tab
  mut_cnv_amp_tab <- mut_cnv_tab[mut_cnv_tab$p_amp > 0 & mut_cnv_tab$p_amp < sig_thres,]
  mut_cnv_amp_tab <- mut_cnv_amp_tab[((mut_cnv_amp_tab$meddiff_bottom_amp > 0) & (mut_cnv_amp_tab$meddiff_upper_amp > 0 ) & mut_cnv_amp_tab$enzyme_type == "kinase") | ((mut_cnv_amp_tab$meddiff_bottom_amp < 0) & (mut_cnv_amp_tab$meddiff_upper_amp < 0 ) & mut_cnv_amp_tab$enzyme_type == "phosphatase"),]
  mut_cnv_amp_tab <- mut_cnv_amp_tab[order(mut_cnv_amp_tab$p_amp),]
  count_enzyme_amp <- data.frame(table(mut_cnv_amp_tab$GENE))
  count_enzyme_amp$genoalt_type <- "amplification"
  
  ## filter amplification-impacted kinase-substrate pairs
  mut_cnv_del_tab <- mut_cnv_tab
  mut_cnv_del_tab <- mut_cnv_tab[mut_cnv_tab$p_del > 0 & mut_cnv_tab$p_del < sig_thres,]
  mut_cnv_del_tab <- mut_cnv_del_tab[((mut_cnv_del_tab$meddiff_bottom_del > 0) & (mut_cnv_del_tab$meddiff_upper_del > 0 ) & mut_cnv_del_tab$enzyme_type == "phosphatase") | ((mut_cnv_del_tab$meddiff_bottom_del < 0) & (mut_cnv_del_tab$meddiff_upper_del < 0 ) & mut_cnv_del_tab$enzyme_type == "kinase"),]
  mut_cnv_del_tab <- mut_cnv_del_tab[order(mut_cnv_del_tab$p_del),]
  count_enzyme_del <- data.frame(table(mut_cnv_del_tab$GENE))
  count_enzyme_del$genoalt_type <- "deletion"
  
  count_enzyme_can <- rbind(count_enzyme_mut, count_enzyme_amp, count_enzyme_del)
  colnames(count_enzyme_can) <- c("GENE", "Freq", "genoalt_type")
  count_enzyme_can <- merge(count_enzyme_can, unique(ptms_site_pairs_sup[, c("GENE", "enzyme_type")]), all.x = T)
  count_enzyme_can$cancer <- cancer
  count_enzyme_cans <- rbind(count_enzyme_cans, count_enzyme_can)
}
count_enzyme_cans <- data.frame(count_enzyme_cans)
count_enzyme_cans$driver_type <- ""
count_enzyme_cans$driver_type[count_enzyme_cans$GENE %in% oncogenes] <- "oncogene"
count_enzyme_cans$driver_type[count_enzyme_cans$GENE %in% tsgs] <- "tsg"

# amp & del - cancer_type_specific ----------------------------------------
for (genoalt_type in c("amplification", "deletion")) {
  tab2p <- count_enzyme_cans
  tab2p <- tab2p[tab2p$genoalt_type == genoalt_type,]
  tab2p <- tab2p[tab2p$Freq > 1,]
  ## divide genes into 2 types: only in one cancer type or >1 cancer types
  count_can_gene <- data.frame(table(tab2p$GENE))
  tab2p$cancer_type_specific <- ""
  tab2p$cancer_type_specific[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq == 1]] <- as.vector(tab2p$cancer[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq == 1]])
  tab2p$cancer_type_specific[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq > 1]] <- "multi-cancer"
  tab2p$cancer_type_specific <- factor(tab2p$cancer_type_specific, levels = c(cancers_sort, "multi-cancer"))
  tab2p$x <- paste0(tab2p$GENE, ":", tab2p$cancer)
  ## sort gene by BRCA on the left, OV in the middle and CO on the right
  
  ## plot cancer type specific ones
  ## highlight oncogene in amplification, tsg in deletion
  if (genoalt_type == "amplification") {
    tab2p$driver_color <- ifelse(tab2p$driver_type == "oncogene", set1[1], "black")
  }
  if (genoalt_type == "deletion") {
    tab2p$driver_color <- ifelse(tab2p$driver_type == "tsg", set1[2], "black")
  }
  tab2p$enzyme_face <- ifelse(tab2p$enzyme_type == "phosphatase", "bold", "plain")
  
  tab2p_sorted <- NULL
  for (cancer in cancers_sort) {
    tab2p_tmp <- tab2p[tab2p$cancer_type_specific == cancer,]
    tab2p_tmp <- tab2p_tmp[order(tab2p_tmp$Freq, decreasing = T),]
    tab2p_tmp <- tab2p_tmp[1:min(top_num, nrow(tab2p_tmp)),]
    tab2p_sorted <- rbind(tab2p_sorted, tab2p_tmp)
  }
  tab2p_sorted$x <- factor(tab2p_sorted$x, levels = as.vector(tab2p_sorted$x))
  
  p <- ggplot(tab2p_sorted)
  p <- p + geom_bar(data=tab2p_sorted, aes(x = x, y = Freq, fill = cancer), stat="identity", position='dodge', color = "#000000")
  p <- p + geom_text(data=tab2p_sorted[tab2p_sorted$Freq > 10,], aes(x = x, y = Freq, color = cancer, label = Freq), nudge_y = 0.05, size = 2)
  # p <- p + facet_grid(.~cancer_type_specific, space = "free", scales = "free")
  p <- p + scale_x_discrete(breaks=as.vector(tab2p_sorted$x), labels=as.vector(tab2p_sorted$GENE))
  p <- p + scale_fill_manual(values = color_cancers2)
  p <- p + scale_color_manual(values = color_cancers2)
  p <- p + scale_y_log10()
  p <- p + theme_nogrid()
  p <- p + theme(axis.title.y = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 0.5, vjust = 0.5, colour = as.vector(tab2p_sorted$driver_color), face = as.vector(tab2p_sorted$enzyme_face)))
  p <- p + theme(panel.spacing.x = unit(units = "line", x = 0))
  p
  fn = paste(makeOutDir(resultD = resultD), 'enzyme_', genoalt_type, '_affect_phosphosites_cancer_type_specific_p_value', sig_thres,'.pdf',sep ="")
  ggsave(filename = fn, width = 6, height = 3)
}

# amp & del - cancer_type_common ----------------------------------------
for (genoalt_type in c("amplification", "deletion")) {
  tab2p <- count_enzyme_cans
  tab2p <- tab2p[tab2p$genoalt_type == genoalt_type,]
  tab2p <- tab2p[tab2p$Freq > 1,]
  ## divide genes into 2 types: only in one cancer type or >1 cancer types
  count_can_gene <- data.frame(table(tab2p$GENE))
  tab2p$cancer_type_specific <- ""
  tab2p$cancer_type_specific[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq == 1]] <- as.vector(tab2p$cancer[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq == 1]])
  tab2p$cancer_type_specific[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq > 1]] <- "multi-cancer"
  tab2p$x <- paste0(tab2p$GENE, ":", tab2p$cancer)
  ## sort gene by BRCA on the left, OV in the middle and CO on the right
  
  ## plot cancer type specific ones
  ## highlight oncogene in amplification, tsg in deletion
  if (genoalt_type == "amplification") {
    tab2p$driver_color <- ifelse(tab2p$driver_type == "oncogene", set1[1], "black")
  }
  if (genoalt_type == "deletion") {
    tab2p$driver_color <- ifelse(tab2p$driver_type == "tsg", set1[2], "black")
  }
  tab2p <- tab2p[tab2p$cancer_type_specific == "multi-cancer",]
  tab2p_sorted <- NULL
  top_genes <- NULL
  for (cancer in cancers_sort) {
    tab2p_tmp <- tab2p[tab2p$cancer_type_specific == "multi-cancer" & tab2p$cancer == cancer,]
    tab2p_tmp <- tab2p_tmp[order(tab2p_tmp$Freq, decreasing = T),]
    top_genes <- c(top_genes, as.vector(tab2p_tmp$GENE[1:min(top_num, nrow(tab2p_tmp))]))
    tab2p_sorted <- rbind(tab2p_sorted, tab2p_tmp)
  }
  top_genes_tab <- data.frame(table(top_genes))
  top_genes_tab <- top_genes_tab[order(top_genes_tab$Freq, decreasing = T),]
  tab2p_sorted <- tab2p_sorted[tab2p_sorted$GENE %in% top_genes_tab$top_genes[1:min(nrow(top_genes_tab), top_num)],]
  tab2p_sorted$cancer <- factor(tab2p_sorted$cancer, levels = cancers_sort)
  tab2p_sorted$GENE <- factor(tab2p_sorted$GENE, levels = as.vector(tab2p_sorted$GENE[!duplicated(tab2p_sorted$GENE)]))
  
  p <- ggplot(tab2p_sorted)
  p <- p + geom_bar(data=tab2p_sorted, aes(x = GENE, y = Freq, fill = cancer), stat="identity", position='stack', color = "#000000")
  p <- p + geom_text(data=tab2p_sorted[tab2p_sorted$Freq > 10,], aes(x = GENE, y = Freq, label = Freq), nudge_y = 0.1, size = 2, color = "#000000")
  p <- p + facet_grid(cancer~., space = "free", scales = "fixed")
  p <- p + scale_fill_manual(values = color_cancers2)
  p <- p + scale_color_manual(values = color_cancers2)
  p <- p + scale_y_log10()
  p <- p + theme_nogrid()
  p <- p + theme(axis.title.y = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 0.5, vjust = 0.5, colour = as.vector(tab2p_sorted$driver_color), face = as.vector(tab2p_sorted$enzyme_face)))
  p <- p + theme(panel.spacing.y = unit(units = "line", x = 0))
  p
  fn = paste(makeOutDir(resultD = resultD), 'enzyme_', genoalt_type, '_affect_phosphosites_cancer_type_shared_p_value', sig_thres,'.pdf',sep ="")
  ggsave(filename = fn, width = 4, height = 3)
}


# barplot for mutations ---------------------------------------------------
for (genoalt_type in c("mutation")) {
  tab2p <- count_enzyme_cans
  tab2p <- tab2p[tab2p$genoalt_type == genoalt_type,]
  tab2p <- tab2p[tab2p$Freq > 1,]
  ## divide genes into 2 types: only in one cancer type or >1 cancer types
  count_can_gene <- data.frame(table(tab2p$GENE))
  tab2p$cancer_type_specific <- ""
  tab2p$cancer_type_specific[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq == 1]] <- as.vector(tab2p$cancer[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq == 1]])
  tab2p$cancer_type_specific[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq > 1]] <- "multi-cancer"
  tab2p$cancer_type_specific <- factor(tab2p$cancer_type_specific, levels = c(cancers_sort, "multi-cancer"))
  tab2p$x <- paste0(tab2p$GENE, ":", tab2p$cancer)
  ## sort gene by BRCA on the left, OV in the middle and CO on the right
  
  ## plot cancer type specific ones
  ## highlight oncogene in amplification, tsg in deletion
  tab2p$driver_color <- "black" 
  tab2p$driver_color[tab2p$driver_type == "oncogene"] <- set1[1]
  tab2p$driver_color[tab2p$driver_type == "tsg"] <- set1[2]
  tab2p$enzyme_face <- ifelse(tab2p$enzyme_type == "phosphatase", "bold", "plain")
  
  tab2p_sorted <- NULL
  for (cancer in cancers_sort) {
    tab2p_tmp <- tab2p[tab2p$cancer == cancer,]
    tab2p_tmp <- tab2p_tmp[order(tab2p_tmp$Freq, decreasing = T),]
    tab2p_sorted <- rbind(tab2p_sorted, tab2p_tmp)
  }
  tab2p_sorted$x <- factor(tab2p_sorted$x, levels = as.vector(tab2p_sorted$x))
  
  p <- ggplot(tab2p_sorted)
  p <- p + geom_bar(data=tab2p_sorted, aes(x = x, y = Freq, fill = cancer), stat="identity", position='dodge', color = "#000000")
  p <- p + geom_text(data=tab2p_sorted[tab2p_sorted$Freq > 10,], aes(x = x, y = Freq,label = Freq), nudge_y = 0.1, size = 2, color = "#000000")
  # p <- p + facet_grid(.~cancer_type_specific, space = "free", scales = "free")
  p <- p + scale_x_discrete(breaks=as.vector(tab2p_sorted$x), labels=as.vector(tab2p_sorted$GENE))
  p <- p + scale_fill_manual(values = color_cancers2)
  p <- p + scale_color_manual(values = color_cancers2)
  p <- p + scale_y_log10()
  p <- p + theme_nogrid()
  p <- p + theme(axis.title.y = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 0.5, vjust = 0.5, colour = as.vector(tab2p_sorted$driver_color), face = as.vector(tab2p_sorted$enzyme_face)))
  p <- p + theme(panel.spacing.x = unit(units = "line", x = 0))
  p
  fn = paste(makeOutDir(resultD = resultD), 'enzyme_', genoalt_type, '_affect_phosphosites_cancer_type_specific_p_value', sig_thres,'.pdf',sep ="")
  ggsave(filename = fn, width = 6, height = 2)
  
  tab2p_driver <- tab2p_sorted[tab2p_sorted$GENE %in% driver_genes$Gene,]
  p <- ggplot(tab2p_driver)
  p <- p + geom_bar(data=tab2p_driver, aes(x = x, y = Freq, fill = cancer), stat="identity", position='dodge', color = "#000000")
  p <- p + geom_text(data=tab2p_driver[tab2p_driver$Freq > 10,], aes(x = x, y = Freq,label = Freq), nudge_y = 0.1, size = 4, color = "#000000")
  # p <- p + facet_grid(.~cancer_type_specific, space = "free", scales = "free")
  p <- p + scale_x_discrete(breaks=as.vector(tab2p_driver$x), labels=as.vector(tab2p_driver$GENE))
  p <- p + scale_fill_manual(values = color_cancers2)
  p <- p + scale_color_manual(values = color_cancers2)
  p <- p + scale_y_log10()
  p <- p + theme_nogrid()
  p <- p + theme(axis.title.y = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 0.5, vjust = 0.5, colour = as.vector(tab2p_driver$driver_color), face = as.vector(tab2p_sorted$enzyme_face)))
  p <- p + theme(panel.spacing.x = unit(units = "line", x = 0))
  p
  fn = paste(makeOutDir(resultD = resultD), 'driver_enzyme_', genoalt_type, '_affect_phosphosites_cancer_type_specific_p_value', sig_thres,'.pdf',sep ="")
  ggsave(filename = fn, width = 4, height = 3)
}


# overview ----------------------------------------------------------------
for (genoalt_type in unique(count_enzyme_cans$genoalt_type)) {
  tab2p <- count_enzyme_cans
  tab2p <- tab2p[tab2p$genoalt_type == genoalt_type,]
  tab2p <- tab2p[tab2p$Freq > 1,]
  ## divide genes into 2 types: only in one cancer type or >1 cancer types
  count_can_gene <- data.frame(table(tab2p$GENE))
  tab2p$cancer_type_specific <- ""
  tab2p$cancer_type_specific[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq == 1]] <- as.vector(tab2p$cancer[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq == 1]])
  tab2p$cancer_type_specific[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq > 1]] <- "multi-cancer"
  tab2p$cancer_type_specific <- factor(tab2p$cancer_type_specific, levels = c(cancers_sort, "multi-cancer"))
  tab2p$x <- paste0(tab2p$GENE, ":", tab2p$cancer)
  tab2p$x <- factor(tab2p$x, levels = as.vector(tab2p$x)[order(tab2p$cancer_type_specific, tab2p$Freq, decreasing = T)])
  
  ## sort gene by BRCA on the left, OV in the middle and CO on the right
  p <- ggplot()
  p <- p + geom_bar(data=tab2p, aes(x = x, y = Freq, fill = cancer), stat="identity", position='dodge', color = "#000000")
  p <- p + geom_text(data=tab2p, aes(x = x, y = Freq, color = cancer, label = Freq), nudge_y = 0.2)
  p <- p + facet_grid(.~cancer_type_specific, space = "free", scales = "free")
  # p <- p + scale_fill_manual(values = c("down" = "#1F78B4", "up" = "#E31A1C"))
  p <- p + scale_x_discrete(breaks=as.vector(tab2p$x), labels=as.vector(tab2p$GENE))
  p <- p + scale_fill_manual(values = color_cancers2)
  p <- p + scale_color_manual(values = color_cancers2)
  p <- p + scale_y_log10()
  p <- p + theme_nogrid()
  p <- p + theme(axis.title.y = element_blank(), axis.ticks.x = element_blank())
  p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
  p <- p + theme(strip.switch.pad.grid = unit(units = "line", x = 0))
  p
  fn = paste(makeOutDir(resultD = resultD), 'enzyme_', genoalt_type, '_affect_phosphosites_p_value', sig_thres,'.pdf',sep ="")
  ggsave(filename = fn, width = 8, height = 3)
}



