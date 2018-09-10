# Yige Wu @ WashU 2018 July
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
sig_thres2test <- c(0.05)
top_num <- 15
networkin_thres <- 5
cancers2process <- cancers_sort
pair_cat2process <- c("shared3can", "uniq_BRCA", "uniq_OV", "uniq_CO")
color_pair <- c("#6A3D9A", color_cancers2[cancers2process])
names(color_pair) <- pair_cat2process
freq_thres <- list(trans = list(high_level = list(amplification = 0, deletion = 0),
                                low_level = list(amplification = 0, deletion = 0)),
                   cis = list(high_level = list(amplification = 0, deletion = 0),
                              low_level = list(amplification = 0, deletion = 0)))
id_vars <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "cancer")
prefix_vars <- c("p", "meddiff_bottom", "meddiff_upper", "meddiff")
color_cna_type <- c(set1[1], set1[2], "#FB9A99", "#A6CEE3"); names(color_cna_type) <- c("high_level_amplification", "high_level_deletion", "low_level_amplification", "low_level_deletion")

# inputs ------------------------------------------------------------------
## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= networkin_thres),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)

## input drive gene list
driver_genes <- read_excel("./Ding_Lab/Projects_Current/TCGA_data/gene_lists/mmc1.xlsx", 
                           sheet = "Table S1", skip = 3)
driver_genes <- data.frame(driver_genes)
oncogenes <- driver_genes$Gene[grepl(x = driver_genes$Tumor.suppressor.or.oncogene.prediction..by.20.20.., pattern = "oncogene")]
tsgs <- driver_genes$Gene[grepl(x = driver_genes$Tumor.suppressor.or.oncogene.prediction..by.20.20.., pattern = "tsg")]
oncogenes <- c(oncogenes, "CTNND1", "PIK3R1")
tsgs <- tsgs[!(tsgs %in% c("CDH1", "CTNND1", "PIK3R1"))]
oncogenes <- unique(oncogenes)
tsgs <- tsgs(oncogenes)

# count per enzyme  ------------------------
for (sig_thres in sig_thres2test) {
  sig_thres <- 0.05
  count_enzyme_cans <- NULL
  mut_cnv_sig_cans <- NULL
  mut_cnv_cans <- NULL

  for (cancer in cancers_sort) {
    mut_cnv_tab <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_cna_impact/", cancer, "_mut_cnv_tab.txt"), data.table = F)
    mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE)
    mut_cnv_tab$pair <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE, ":", mut_cnv_tab$SUB_MOD_RSD)
    mut_cnv_tab <- merge(mut_cnv_tab, unique(ptms_site_pairs_sup[, c("GENE", "enzyme_type")]), all.x = T)
    
    ##filter out bad annotated enzyme type
    mut_cnv_tab <- mut_cnv_tab[!(mut_cnv_tab$GENE == "PTPRG" & mut_cnv_tab$enzyme_type == "kinase"),]
    
    ## annotate cis and trans
    mut_cnv_tab$SELF <- "trans"
    mut_cnv_tab$SELF[as.vector(mut_cnv_tab$GENE) == as.vector(mut_cnv_tab$SUB_GENE)] <- "cis"
    
    ## filter mutation-impacted kinase-substrate pairs
    mut_tab <- mut_cnv_tab
    mut_tab <- mut_tab[mut_tab$p_mut > 0,]
    
    sig_mut_tab <- mut_tab[mut_tab$p_mut < sig_thres,]
    sig_mut_tab <- sig_mut_tab[(sig_mut_tab$meddiff_bottom_mut > 0) == (sig_mut_tab$meddiff_upper_mut > 0),]
    sig_mut_tab <- sig_mut_tab[order(sig_mut_tab$p_mut),]
    
    shallow_amp_tab <- mut_cnv_tab
    shallow_amp_tab <- shallow_amp_tab[shallow_amp_tab$p_shallow_amp > 0,]
    sig_shallow_amp_tab <- shallow_amp_tab[shallow_amp_tab$p_shallow_amp < sig_thres,]
    sig_shallow_amp_tab <- sig_shallow_amp_tab[((sig_shallow_amp_tab$meddiff_bottom_shallow_amp > 0) & (sig_shallow_amp_tab$meddiff_upper_shallow_amp > 0 ) & sig_shallow_amp_tab$enzyme_type == "kinase") | ((sig_shallow_amp_tab$meddiff_bottom_shallow_amp < 0) & (sig_shallow_amp_tab$meddiff_upper_shallow_amp < 0 ) & sig_shallow_amp_tab$enzyme_type == "phosphatase"),]
    sig_shallow_amp_tab <- sig_shallow_amp_tab[order(sig_shallow_amp_tab$p_shallow_amp),]
    
    deep_amp_tab <- mut_cnv_tab
    deep_amp_tab <- deep_amp_tab[deep_amp_tab$p_deep_amp > 0,]
    sig_deep_amp_tab <- deep_amp_tab[deep_amp_tab$p_deep_amp < sig_thres,]
    sig_deep_amp_tab <- sig_deep_amp_tab[((sig_deep_amp_tab$meddiff_bottom_deep_amp > 0) & (sig_deep_amp_tab$meddiff_upper_deep_amp > 0 ) & sig_deep_amp_tab$enzyme_type == "kinase") | ((sig_deep_amp_tab$meddiff_bottom_deep_amp < 0) & (sig_deep_amp_tab$meddiff_upper_deep_amp < 0 ) & sig_deep_amp_tab$enzyme_type == "phosphatase"),]
    sig_deep_amp_tab <- sig_deep_amp_tab[order(sig_deep_amp_tab$p_deep_amp),]
    
    ## filter deletion-impacted kinase-substrate pairs
    shallow_del_tab <- mut_cnv_tab
    shallow_del_tab <- shallow_del_tab[shallow_del_tab$p_shallow_del > 0,]
    sig_shallow_del_tab <- shallow_del_tab[shallow_del_tab$p_shallow_del < sig_thres,]
    sig_shallow_del_tab <- sig_shallow_del_tab[((sig_shallow_del_tab$meddiff_bottom_shallow_del > 0) & (sig_shallow_del_tab$meddiff_upper_shallow_del > 0 ) & sig_shallow_del_tab$enzyme_type == "phosphatase") | ((sig_shallow_del_tab$meddiff_bottom_shallow_del < 0) & (sig_shallow_del_tab$meddiff_upper_shallow_del < 0 ) & sig_shallow_del_tab$enzyme_type == "kinase"),]
    
    ## filter amplification-impacted kinase-substrate pairs
    deep_del_tab <- mut_cnv_tab
    deep_del_tab <- deep_del_tab[deep_del_tab$p_deep_del > 0,]
    sig_deep_del_tab <- deep_del_tab[deep_del_tab$p_deep_del < sig_thres,]
    sig_deep_del_tab <- sig_deep_del_tab[((sig_deep_del_tab$meddiff_bottom_deep_del > 0) & (sig_deep_del_tab$meddiff_upper_deep_del > 0 ) & sig_deep_del_tab$enzyme_type == "phosphatase") | ((sig_deep_del_tab$meddiff_bottom_deep_del < 0) & (sig_deep_del_tab$meddiff_upper_deep_del < 0 ) & sig_deep_del_tab$enzyme_type == "kinase"),]
    
    ## count the occurence of each enzyme
    count_enzyme_mut <- data.frame(table(sig_mut_tab[, c("GENE", "SELF")]))
    if (nrow(count_enzyme_mut) > 0){
      count_enzyme_mut$genoalt_type <- "mutation"
      
      ## format the mutation table
      mut_tab.f <- mut_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "mut"), "pair")]
      colnames(mut_tab.f) <- c(id_vars, prefix_vars, "pair")
      mut_tab.f$genoalt_type <- "mutation"
      sig_mut_tab.f <- mut_tab.f[mut_tab.f$pair %in% sig_mut_tab$pair,]
    }
    
    count_enzyme_shallow_amp <- data.frame(table(sig_shallow_amp_tab[, c("GENE", "SELF")]))
    if (nrow(count_enzyme_shallow_amp) > 0){
      count_enzyme_shallow_amp$genoalt_type <- "low_level_amplification"
      shallow_amp_tab.f <- shallow_amp_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "shallow_amp"), "pair")]
      colnames(shallow_amp_tab.f) <- c(id_vars, prefix_vars, "pair")
      shallow_amp_tab.f$genoalt_type <- "low_level_amplification"
      
      sig_shallow_amp_tab.f <- shallow_amp_tab.f[shallow_amp_tab.f$pair %in% sig_shallow_amp_tab$pair, ]
    }
    
    
    count_enzyme_shallow_del <- data.frame(table(sig_shallow_del_tab[, c("GENE", "SELF")]))
    if (nrow(count_enzyme_shallow_del) > 0){
      count_enzyme_shallow_del$genoalt_type <- "low_level_deletion"
      shallow_del_tab.f<- shallow_del_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "shallow_del"), "pair")]
      colnames(shallow_del_tab.f) <- c(id_vars, prefix_vars, "pair")
      shallow_del_tab.f$genoalt_type <- "low_level_deletion"
      
      sig_shallow_del_tab.f <- shallow_del_tab.f[shallow_del_tab.f$pair %in% sig_shallow_del_tab$pair,]
    }
    count_enzyme_deep_amp <- data.frame(table(sig_deep_amp_tab[, c("GENE", "SELF")]))
    if (nrow(count_enzyme_deep_amp) > 0){
      count_enzyme_deep_amp$genoalt_type <- "high_level_amplification"
      deep_amp_tab.f<- deep_amp_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "deep_amp"), "pair")]
      colnames(deep_amp_tab.f) <- c(id_vars, prefix_vars, "pair")
      deep_amp_tab.f$genoalt_type <- "high_level_amplification"
      
      sig_deep_amp_tab.f <- deep_amp_tab.f[deep_amp_tab.f$pair %in% sig_deep_amp_tab$pair,]
      
    }
    
    count_enzyme_deep_del <- data.frame(table(sig_deep_del_tab[, c("GENE", "SELF")]))
    if (nrow(count_enzyme_deep_del) > 0){
      count_enzyme_deep_del$genoalt_type <- "high_level_deletion"
      deep_del_tab.f<- deep_del_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "deep_del"), "pair")]
      colnames(deep_del_tab.f) <- c(id_vars, prefix_vars, "pair")
      deep_del_tab.f$genoalt_type <- "high_level_deletion"
      
      sig_deep_del_tab.f <- deep_del_tab.f[deep_del_tab.f$pair %in% sig_deep_del_tab$pair,]
    }
    
    
    count_enzyme_can <- rbind(count_enzyme_mut, count_enzyme_shallow_amp, count_enzyme_shallow_del, count_enzyme_deep_amp, count_enzyme_deep_del)
    colnames(count_enzyme_can) <- c("GENE", "SELF", "Freq", "genoalt_type")
    count_enzyme_can <- merge(count_enzyme_can, unique(ptms_site_pairs_sup[, c("GENE", "enzyme_type")]), all.x = T)
    count_enzyme_can$cancer <- cancer
    count_enzyme_cans <- rbind(count_enzyme_cans, count_enzyme_can)
    
    mut_cnv_sig_can <- rbind(sig_mut_tab.f, sig_deep_del_tab.f, sig_deep_amp_tab.f, sig_shallow_del_tab.f, sig_shallow_amp_tab.f)
    mut_cnv_sig_cans <- rbind(mut_cnv_sig_cans, mut_cnv_sig_can)
    mut_cnv_sig_cans <- unique(mut_cnv_sig_cans)
    
    mut_cnv_can <- rbind(mut_tab.f, deep_del_tab.f, deep_amp_tab.f, shallow_del_tab.f, shallow_amp_tab.f)
    mut_cnv_cans <- rbind(mut_cnv_cans, mut_cnv_can)
    mut_cnv_cans <- unique(mut_cnv_cans)
  }
  count_enzyme_cans <- data.frame(count_enzyme_cans)
  count_enzyme_cans$driver_type <- ""
  count_enzyme_cans$driver_type[count_enzyme_cans$GENE %in% oncogenes] <- "oncogene"
  count_enzyme_cans$driver_type[count_enzyme_cans$GENE %in% tsgs] <- "tsg"
  
  mut_cnv_sig_cans$SELF <- ifelse(as.vector(mut_cnv_sig_cans$GENE) == as.vector(mut_cnv_sig_cans$SUB_GENE), "cis", "trans")
  mut_cnv_cans$SELF <- ifelse(as.vector(mut_cnv_cans$GENE) == as.vector(mut_cnv_cans$SUB_GENE), "cis", "trans")

  write.table(x = count_enzyme_cans, file = paste0(makeOutDir(resultD = resultD), "count_enzyme_cans_sig_thres", sig_thres, ".txt"), row.names = F, quote = F, sep = "\t")
  write.table(x = mut_cnv_sig_cans, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_sig_cans_sig_thres", sig_thres, ".txt"), row.names = F, quote = F, sep = "\t")
  write.table(x = mut_cnv_cans, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_sig_cans.txt"), row.names = F, quote = F, sep = "\t")

  for (self in c("trans", "cis")) {
    count_enzyme_cans_self <- count_enzyme_cans[count_enzyme_cans$SELF == self,]

    # barplot for mutations ---------------------------------------------------
    # for (genoalt_type in c("mutation")) {
    #   tab2p <- count_enzyme_cans_self
    #   tab2p <- tab2p[tab2p$genoalt_type == genoalt_type,]
    #   tab2p <- tab2p[tab2p$Freq > 0,]
    #   ## divide genes into 2 types: only in one cancer type or >1 cancer types
    #   count_can_gene <- data.frame(table(tab2p$GENE))
    #   tab2p$cancer_type_specific <- ""
    #   tab2p$cancer_type_specific[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq == 1]] <- as.vector(tab2p$cancer[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq == 1]])
    #   tab2p$cancer_type_specific[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq > 1]] <- "shared3can"
    #   tab2p$cancer_type_specific <- factor(tab2p$cancer_type_specific, levels = c(cancers_sort, "shared3can"))
    #   tab2p$x <- paste0(tab2p$GENE, ":", tab2p$cancer)
    #   ## sort gene by BRCA on the left, OV in the middle and CO on the right
    #   
    #   ## plot cancer type specific ones
    #   ## highlight oncogene in amplification, tsg in deletion
    #   tab2p$driver_color <- "black" 
    #   tab2p$driver_color[tab2p$driver_type == "oncogene"] <- set1[1]
    #   tab2p$driver_color[tab2p$driver_type == "tsg"] <- set1[2]
    #   tab2p$enzyme_face <- ifelse(tab2p$enzyme_type == "phosphatase", "bold", "plain")
    #   
    #   tab2p_sorted <- NULL
    #   for (cancer in cancers_sort) {
    #     tab2p_tmp <- tab2p[tab2p$cancer == cancer,]
    #     tab2p_tmp <- tab2p_tmp[order(tab2p_tmp$Freq, decreasing = T),]
    #     tab2p_sorted <- rbind(tab2p_sorted, tab2p_tmp)
    #   }
    #   tab2p_sorted$x <- factor(tab2p_sorted$x, levels = as.vector(tab2p_sorted$x))
    #   
    #   p <- ggplot(tab2p_sorted)
    #   p <- p + geom_bar(data=tab2p_sorted, aes(x = x, y = Freq, fill = cancer), stat="identity", position='dodge', color = "#000000")
    #   if (self == "cis") {
    #     p <- p + geom_text(data=tab2p_sorted[tab2p_sorted$Freq > 0,], aes(x = x, y = Freq,label = Freq), nudge_y = 0.1, size = 2, color = "#000000")
    #   } else {
    #     p <- p + geom_text(data=tab2p_sorted[tab2p_sorted$Freq > 10,], aes(x = x, y = Freq,label = Freq), nudge_y = 0.1, size = 2, color = "#000000")
    #   }
    #   p <- p + scale_x_discrete(breaks=as.vector(tab2p_sorted$x), labels=as.vector(tab2p_sorted$GENE))
    #   p <- p + scale_fill_manual(values = color_cancers2)
    #   p <- p + scale_color_manual(values = color_cancers2)
    #   p <- p + scale_y_log10()
    #   p <- p + theme_nogrid()
    #   p <- p + theme(axis.title.y = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
    #   p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
    #   p <- p + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 0.5, vjust = 0.5, colour = as.vector(tab2p_sorted$driver_color), face = as.vector(tab2p_sorted$enzyme_face)))
    #   p <- p + theme(panel.spacing.x = unit(units = "line", x = 0))
    #   p
    #   fn = paste(makeOutDir(resultD = resultD), 'enzyme_', genoalt_type, '_affect_', self, '_phosphosites_cancer_type_specific_p_value', sig_thres,'.pdf',sep ="")
    #   ggsave(filename = fn, width = 6, height = 2)
    #   
    #   tab2p_driver <- tab2p_sorted[tab2p_sorted$GENE %in% driver_genes$Gene,]
    #   p <- ggplot(tab2p_driver)
    #   p <- p + geom_bar(data=tab2p_driver, aes(x = x, y = Freq, fill = cancer), stat="identity", position='dodge', color = "#000000")
    #   p <- p + geom_text(data=tab2p_driver[tab2p_driver$Freq > 0,], aes(x = x, y = Freq,label = Freq), nudge_y = 0.1, size = 4, color = "#000000")
    #   # p <- p + facet_grid(.~cancer_type_specific, space = "free", scales = "free")
    #   p <- p + scale_x_discrete(breaks=as.vector(tab2p_driver$x), labels=as.vector(tab2p_driver$GENE))
    #   p <- p + scale_fill_manual(values = color_cancers2)
    #   p <- p + scale_color_manual(values = color_cancers2)
    #   p <- p + scale_y_log10()
    #   p <- p + theme_nogrid()
    #   p <- p + theme(axis.title.y = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
    #   p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
    #   p <- p + theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 0.5, vjust = 0.5, colour = as.vector(tab2p_driver$driver_color), face = as.vector(tab2p_sorted$enzyme_face)))
    #   p <- p + theme(panel.spacing.x = unit(units = "line", x = 0))
    #   p
    #   fn = paste(makeOutDir(resultD = resultD), 'driver_enzyme_', genoalt_type, '_affect_', self, '_phosphosites_cancer_type_specific_p_value', sig_thres,'.pdf',sep ="")
    #   ggsave(filename = fn, width = ifelse(self == "cis", 2, 4), height = 3)
    #   
    #   tab2p_driver <- tab2p_sorted[(tab2p_sorted$GENE %in% driver_genes$Gene & tab2p_sorted$cancer == "CO") | tab2p_sorted$cancer != "CO",]
    #   p <- ggplot(tab2p_driver)
    #   p <- p + geom_bar(data=tab2p_driver, aes(x = x, y = Freq, fill = cancer), stat="identity", position='dodge', color = "#000000")
    #   p <- p + geom_text(data=tab2p_driver[tab2p_driver$Freq > 0,], aes(x = x, y = Freq,label = Freq), nudge_y = 0.1, size = 4, color = "#000000")
    #   # p <- p + facet_grid(.~cancer_type_specific, space = "free", scales = "free")
    #   p <- p + scale_x_discrete(breaks=as.vector(tab2p_driver$x), labels=as.vector(tab2p_driver$GENE))
    #   p <- p + scale_fill_manual(values = color_cancers2)
    #   p <- p + scale_color_manual(values = color_cancers2)
    #   p <- p + scale_y_log10()
    #   p <- p + theme_nogrid()
    #   p <- p + theme(axis.title.y = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
    #   p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
    #   p <- p + theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 0.5, vjust = 0.5, colour = as.vector(tab2p_driver$driver_color), face = as.vector(tab2p_sorted$enzyme_face)))
    #   p <- p + theme(panel.spacing.x = unit(units = "line", x = 0))
    #   p
    #   fn = paste(makeOutDir(resultD = resultD), 'driver_enzyme_', genoalt_type, '_affect_', self, '_phosphosites_cancer_type_specific_p_value', sig_thres,'_V2.pdf',sep ="")
    #   ggsave(filename = fn, width = 5, height = 3)
    # }
    for (genoalt_type in c("amplification", "deletion")) {
      tab2p_sorted <- NULL
      for (cna_level in c("high_level", "low_level")) {
      # cna_level <- "high_level"
        # genoalt_type <- "amplification"
        tab2p <- count_enzyme_cans_self
        tab2p <- tab2p[tab2p$genoalt_type == paste0(cna_level, "_", genoalt_type),]
        tab2p <- tab2p[tab2p$Freq > freq_thres[[self]][[cna_level]][[genoalt_type]],]
        ## divide genes into 2 types: only in one cancer type or >1 cancer types
        if (nrow(tab2p) > 0) {
          count_can_gene <- data.frame(table(tab2p$GENE))
          tab2p$x <- paste0(tab2p$GENE, ":", tab2p$cancer, ":", tab2p$genoalt_type)
          
          ## plot cancer type specific ones
          ## highlight oncogene in amplification, tsg in deletion
          if (genoalt_type == "amplification") {
            tab2p$driver_color <- ifelse(tab2p$driver_type == "oncogene", set1[1], "black")
          }
          if (genoalt_type == "deletion") {
            tab2p$driver_color <- ifelse(tab2p$driver_type == "tsg", set1[2], "black")
          }
          tab2p$enzyme_face <- ifelse(tab2p$enzyme_type == "phosphatase", "bold", "plain")
          # amp & del - cancer_type_specific ----------------------------------------
          for (cancer in cancers_sort) {
            tab2p_tmp <- tab2p[tab2p$cancer == cancer & tab2p$GENE %in% CNAs[[cancer]][[genoalt_type]],]
            if (nrow(tab2p_tmp) > 0){
              tab2p_tmp <- tab2p_tmp[order(tab2p_tmp$Freq, decreasing = T),]
              tab2p_tmp <- tab2p_tmp[1:min(top_num, nrow(tab2p_tmp)),]
              tab2p_sorted <- rbind(tab2p_sorted, tab2p_tmp)
            }
          }
        }
      }
      tab2p_sorted$x <- factor(tab2p_sorted$x, levels = as.vector(tab2p_sorted$x))
      tab2p_sorted$cancer <- factor(tab2p_sorted$cancer, levels = cancers_sort)
      p <- ggplot(tab2p_sorted)
      p <- p + geom_bar(data=tab2p_sorted, aes(x = GENE, y = Freq, fill = genoalt_type), stat="identity", position='dodge', color = "#000000")
      p <- p + facet_grid(genoalt_type~cancer, space = "free", scales = "free_x")
      p <- p + theme_nogrid()
      if (self == "cis") {
        if (genoalt_type == "deletion") {
          p <- p + geom_text(data=tab2p_sorted[tab2p_sorted$Freq > 0,], aes(x = GENE, y = Freq,label = Freq), nudge_y = 0.3, size = 3, color = "#000000")
        } else {
          p <- p + geom_text(data=tab2p_sorted[tab2p_sorted$Freq > 0,], aes(x = GENE, y = Freq,label = Freq), nudge_y = 1, size = 2, color = "#000000")
        }
      } else {
        p <- p + geom_text(data=tab2p_sorted[tab2p_sorted$Freq > 0,], aes(x = GENE, y = Freq, label = Freq), nudge_y = 0.1, size = 2, color = "#000000")
      }
      # p <- p + scale_x_discrete(breaks=as.vector(tab2p_sorted$x), labels=as.vector(tab2p_sorted$GENE))
      p <- p + scale_fill_manual(values = color_cna_type)
      if (self == "trans") {
        p <- p + scale_y_log10()
      }
      p <- p + theme(strip.text.y = element_text(size = ifelse(self == "cis", 5, 10)))
      p <- p + theme(axis.title.y = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
      p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
      p <- p + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 0.5, vjust = 0.5))
      p <- p + theme(panel.spacing.x = unit(units = "line", x = 0))
      p
      fn = paste(makeOutDir(resultD = resultD), 'enzyme_', genoalt_type, '_affect_', self, '_phosphosites_cancer_type_specific_CNAs_p_value', sig_thres,'.pdf',sep ="")
      ggsave(filename = fn, width = ifelse(self == "cis", 4, 6), height =  ifelse(self == "cis", 3, 4))
    }
  }
}
