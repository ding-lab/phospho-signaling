# Yige Wu @ WashU 2018 Mar
# annotate table f hot and cold kinase-substrate pairs by fisher's test

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(dplyr)

# set variables -----------------------------------------------------------
var_list <- c('pro_kin', 'pho_kin')
names(var_list) <- c('cis', 'trans')
subtypes_list <- c('pam50', 'MSI_type', 'tumor_normal_WBasal'); names(subtypes_list) <- c("BRCA", "CO", "OV")
subtype_sort = c("BRCA_LumA", "BRCA_LumB", "BRCA_Her2", "BRCA_Basal" ,"BRCA_Normal", "OV_Tumor","OV_Normal","CO_MSI","CO_MSS","CO_other","CO_Normal")
myvjust = 0.7
myhjust = 0
cancers_sort <- c("BRCA", "OV", "CO")
diffexp_sig <- 0.2

# inputs ------------------------------------------------------------------
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
sup_cans_tab_pp <- sup_cans_tab[sup_cans_tab$enzyme_type == enzyme_type,]

kin_acc_ids <- fread(input = "./cptac2p/analysis_results/enzyme_acc_ids.txt", data.table = F)
# load drivers for all cancers -------------------------------------------------
driver_cans <- NULL
for (cancer in cancers_sort) {
  tmp <- loadGeneList(gene_type = "driver", cancer = cancer, is.soft.limit = "soft")
  driver_cans <- c(driver_cans, tmp)
}

driver_cans <- unique(driver_cans)


# load cancer hallmark pathways -------------------------------------------
load("~/Box Sync/pan3can_shared_data/Gene_family/2015-08-01_Gene_Set.RData")
keywords <- c("cell cycle", "PI3K","Akt","MAPK","adherens","Mismatch","EMT", "apoptosis","immunological","stromal","transmembrane","receptors","integrin",
              "TGFβ","LKB1","AMPK","TSC","mTOR","Ras","Notch","Wnt", "catenin",  "p53","RTK","erbb", "Spliceosome")
all_keggs <- names(KEGG)
keggs <- NULL
for (key in keywords) {
  keggs <- c(keggs, all_keggs[grepl(key, all_keggs, ignore.case = T)])
}
keggs <- unique(keggs)

for (diffexp_dir in c("up")) {
  self <- 'trans'
  
  ## input kinases functional annotation
  func_anno_tab <- fread(input = paste0(ppnD, "diffexp/tables/high&low_sub_in_high&low_es_pairs/", diffexp_dir, "/", diffexp_dir, "_kinase_DAVID.txt"), data.table = F)
  func_anno_tab <- merge(func_anno_tab, kin_acc_ids, by.x = c("ID"), by.y = c("KIN_ACC_ID"), all.x = T)
  rownames(func_anno_tab) <- func_anno_tab$GENE
  
  t0 <- sup_cans_tab_pk
  t0 <- t0[!is.na(t0$KSEA_pvalue) & !is.na(t0$substrate_log2FC),]
  t0 <- t0[!is.na(t0$substrate_direction),]
  t0 <- t0[t0$substrate_direction == diffexp_dir,]
  t0 <- t0[abs(as.vector(t0$substrate_log2FC)) > 1,]
  t0 <- t0[t0$KSEA_pvalue < diffexp_sig,]
  t0 <- t0[t0$KSEA_enzyme_direction == diffexp_dir,]
  enzymes <- unique(t0$GENE)
  
  var <- var_list[self]
  fdr_var <- paste("FDR_",var,sep = "");  coef_var <- paste("coef_",var,sep = "")
  
  if (nrow(t0) > 0){
    t1 <- t0
    
    ## draw pairs with given kinase/substrate list
    if (TRUE){
      t2 <- unique(t1[,c("GENE", 'KSEA_log2FC', "Cancer")])
      
      lim = max(abs(max(t2$KSEA_log2FC, na.rm = T)),abs(min(t2$KSEA_log2FC, na.rm = T)))
      cap <- min(2, lim)
      t2$KSEA_log2FC_capped <- t2$KSEA_log2FC
      t2$KSEA_log2FC_capped[t2$KSEA_log2FC > cap] <- cap
      t2$KSEA_log2FC_capped[t2$KSEA_log2FC < (-cap)] <- (-cap)
      
      ## get the number of consistently high pairs
      cons_tab <- data.frame(table(unique(t0[, c("GENE", "pair", "Cancer")])[, c("GENE", "Cancer")]))
      t2 <- merge(t2, cons_tab, all.x = T)
      
      ## sort cancers
      t2$Cancer <- factor(t2$Cancer, levels = cancers_sort)
      ## sort kinases
      cons_tab_cans <- data.frame(table(unique(t0[, c("GENE", "pair", "Cancer")])[, c("GENE")]))
      colnames(cons_tab_cans) <- c("GENE", "Freq_cans")
      t2 <- merge(t2, cons_tab_cans, all.x = T)
      t2 <- t2[t2$GENE %in% t2$GENE[t2$Freq >= 4],]
      
      ## annotate kinases to pathways
      enzyme_tab <- data.frame(enzyme = unique(t1[,c("GENE")]))
      for (kegg in keggs) {
        kegg_id <- str_split(string = kegg, pattern = "\\t")[[1]][1]
        kegg_col <- func_anno_tab[as.vector(enzyme_tab$enzyme), "KEGG_PATHWAY"]
        enzyme_tab[, kegg] <- grepl(pattern = kegg_id, x = kegg_col, ignore.case = T)
      }
      enzyme_tab$hallmark_count <- rowSums(enzyme_tab[, keggs])
      hallmark_assign <- vector(mode = "character", length = nrow(enzyme_tab))
      for (i in 1:nrow(enzyme_tab)) {
        if (enzyme_tab$hallmark_count[i] == 0) {
          hallmark_assign[i] <- "other"
        }
        if (enzyme_tab$hallmark_count[i] > 0) {
          hallmark_assign[i] <- keggs[min(which(enzyme_tab[i, keggs] == T))]
        }
      }
      enzyme_tab$hallmark_assign <- hallmark_assign
      
      t2 <- merge(t2, enzyme_tab[, c("enzyme", "hallmark_assign")], by.x = c('GENE'), by.y = c('enzyme'), all.x = T)
      ## order kinases by their number of consistent high/low pairs
      t2$GENE <- factor(t2$GENE, levels = as.vector(cons_tab_cans$GENE)[order(cons_tab_cans$Freq_cans, decreasing = F)])
      t2$hallmark_assign <- factor(t2$hallmark_assign, levels = c(keggs, "other"))
      plot3 = ggplot(t2)
      plot3 = plot3 + geom_tile(aes(x=Cancer, y=GENE, fill=KSEA_log2FC_capped), color=NA)#, linetype="blank")
      plot3 = plot3 + scale_fill_gradientn(name= "KSEA_log2FC", na.value=NA, colours=RdBu1024, limits = c(-lim, lim))
      plot3 = plot3 + geom_text(data = t2[t2$Freq > 0,], mapping = aes(x=Cancer, y=GENE, label = Freq))
      plot3 = plot3 + facet_grid(hallmark_assign~., drop=T, space = "free",scales = "free")
      plot3 = plot3 + theme_bw() + theme_nogrid()
      plot3 = plot3 + theme(axis.title = element_blank(), axis.text.y = element_text(size = 10), strip.text.y = element_text(size = 5))
      plot3 = plot3 + theme(axis.text.x = element_text(angle = -30, vjust = myvjust, hjust = myhjust))
      plot3 = plot3 + theme(axis.ticks=element_blank(), legend.position="bottom")
      plot3 = plot3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
      fn = paste(makeOutDir(resultD = resultD), diffexp_dir, '_', enzyme_type ,'substrate_pairs_number.pdf',sep ="")
      plot3
      ggsave(filename = fn, width = 2, height = 12*nrow(t2)/120)
    } else {
      print('not enough data')
    }
  }
}

for (diffexp_dir in c("down")) {
  self <- 'trans'
  
  ## input kinases functional annotation
  func_anno_tab <- fread(input = paste0(ppnD, "diffexp/tables/high&low_sub_in_high&low_es_pairs/", diffexp_dir, "/", diffexp_dir, "_kinase_DAVID.txt"), data.table = F)
  func_anno_tab <- merge(func_anno_tab, kin_acc_ids, by.x = c("ID"), by.y = c("KIN_ACC_ID"), all.x = T)
  rownames(func_anno_tab) <- func_anno_tab$GENE
  
  t0 <- sup_cans_tab_pk
  t0 <- t0[!is.na(t0$KSEA_pvalue) & !is.na(t0$substrate_log2FC),]
  t0 <- t0[!is.na(t0$substrate_direction),]
  t0 <- t0[t0$substrate_direction == diffexp_dir,]
  t0 <- t0[abs(as.vector(t0$substrate_log2FC)) > 1,]
  t0 <- t0[t0$KSEA_pvalue < diffexp_sig,]
  t0 <- t0[t0$KSEA_enzyme_direction == diffexp_dir,]
  enzymes <- unique(t0$GENE)
  
  var <- var_list[self]
  fdr_var <- paste("FDR_",var,sep = "");  coef_var <- paste("coef_",var,sep = "")
  
  if (nrow(t0) > 0){
    t1 <- t0
    
    ## draw pairs with given kinase/substrate list
    if (TRUE){
      t2 <- unique(t1[,c("GENE", 'KSEA_log2FC', "Cancer")])
      
      lim = max(abs(max(t2$KSEA_log2FC, na.rm = T)),abs(min(t2$KSEA_log2FC, na.rm = T)))
      cap <- min(2, lim)
      t2$KSEA_log2FC_capped <- t2$KSEA_log2FC
      t2$KSEA_log2FC_capped[t2$KSEA_log2FC > cap] <- cap
      t2$KSEA_log2FC_capped[t2$KSEA_log2FC < (-cap)] <- (-cap)
      
      ## get the number of consistently high pairs
      cons_tab <- data.frame(table(unique(t0[, c("GENE", "pair", "Cancer")])[, c("GENE", "Cancer")]))
      t2 <- merge(t2, cons_tab, all.x = T)
      
      ## sort cancers
      t2$Cancer <- factor(t2$Cancer, levels = cancers_sort)
      ## sort kinases
      cons_tab_cans <- data.frame(table(unique(t0[, c("GENE", "pair", "Cancer")])[, c("GENE")]))
      colnames(cons_tab_cans) <- c("GENE", "Freq_cans")
      t2 <- merge(t2, cons_tab_cans, all.x = T)
      t2 <- t2[t2$GENE %in% t2$GENE[t2$Freq >= 5],]
      
      ## annotate kinases to pathways
      enzyme_tab <- data.frame(enzyme = unique(t1[,c("GENE")]))
      for (kegg in keggs) {
        kegg_id <- str_split(string = kegg, pattern = "\\t")[[1]][1]
        kegg_col <- func_anno_tab[as.vector(enzyme_tab$enzyme), "KEGG_PATHWAY"]
        enzyme_tab[, kegg] <- grepl(pattern = kegg_id, x = kegg_col, ignore.case = T)
      }
      enzyme_tab$hallmark_count <- rowSums(enzyme_tab[, keggs])
      hallmark_assign <- vector(mode = "character", length = nrow(enzyme_tab))
      for (i in 1:nrow(enzyme_tab)) {
        if (enzyme_tab$hallmark_count[i] == 0) {
          hallmark_assign[i] <- "other"
        }
        if (enzyme_tab$hallmark_count[i] > 0) {
          hallmark_assign[i] <- keggs[min(which(enzyme_tab[i, keggs] == T))]
        }
      }
      enzyme_tab$hallmark_assign <- hallmark_assign
      
      t2 <- merge(t2, enzyme_tab[, c("enzyme", "hallmark_assign")], by.x = c('GENE'), by.y = c('enzyme'), all.x = T)
      t2 <- t2[t2$hallmark_assign != "other",]
      ## order kinases by their number of consistent high/low pairs
      t2$GENE <- factor(t2$GENE, levels = as.vector(cons_tab_cans$GENE)[order(cons_tab_cans$Freq_cans, decreasing = F)])
      t2$hallmark_assign <- factor(t2$hallmark_assign, levels = c(keggs, "other"))
      plot3 = ggplot(t2)
      plot3 = plot3 + geom_tile(aes(x=Cancer, y=GENE, fill=KSEA_log2FC_capped), color=NA)#, linetype="blank")
      plot3 = plot3 + scale_fill_gradientn(name= "KSEA_log2FC", na.value=NA, colours=RdBu1024, limits = c(-lim, lim))
      plot3 = plot3 + geom_text(data = t2[t2$Freq > 0,], mapping = aes(x=Cancer, y=GENE, label = Freq))
      plot3 = plot3 + facet_grid(hallmark_assign~., drop=T, space = "free",scales = "free")
      plot3 = plot3 + theme_bw() + theme_nogrid()
      plot3 = plot3 + theme(axis.title = element_blank(), axis.text.y = element_text(size = 10), strip.text.y = element_text(size = 5))
      plot3 = plot3 + theme(axis.text.x = element_text(angle = -30, vjust = myvjust, hjust = myhjust))
      plot3 = plot3 + theme(axis.ticks=element_blank(), legend.position="bottom")
      plot3 = plot3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
      fn = paste(makeOutDir(resultD = resultD), diffexp_dir, '_', enzyme_type ,'substrate_pairs_number.pdf',sep ="")
      plot3
      ggsave(filename = fn, width = 2, height = 12*nrow(t2)/120)
    } else {
      print('not enough data')
    }
  }
}





