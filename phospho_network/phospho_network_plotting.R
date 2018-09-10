# Yige Wu @ WashU 2018 Jan
## shared plotting functions for phospho_network

source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source(paste0(baseD, "cptac2p_analysis/pan3can_aes.R")) # aes for general purposes; it should be one directory out of the working directory

# library -----------------------------------------------------------------
library(ggplot2)
library(gtable)

# static parameters -------------------------------------------------------
prefix <- NULL
prefix["BRCA"] <- "BRCA_BI"
prefix["OV"] <- "OV_PNNL"
prefix["CO"] <- "CO_PNNL"
sig <- 0.05 # significance level
setwd(baseD)

color_cancers <- c(colors['BRCA'], colors['OV'], colors['COAD'], "#bdbdbd"); names(color_cancers) <- c("BRCA", "OV", "CO", "other")
color_cancers2 <- c(colors['BRCA'], colors['OV'],  "#39BEB1", "#fb9a99", "#bdbdbd"); names(color_cancers2) <- c("BRCA", "OV", "CO", "UCEC","other")

## colors for different amino acid PTMs
color_aa <- c("#8da0cb", "#e78ac3", "#a6d854"); names(color_aa) <- c("S", "T", "Y")

bubble_heatmap = function(cancer_id, protein, self, id, table_bubble, table_substrate, table_kinase, outDir, subtype_sort, h = 12, w = 5.5) {
  ## Usage: plut bubble plot and heatmap side by side
  ## input:
  ##        identifier(string): cancer type of the dataset, kinase/phophotase substrate pairs; cis/trans; additional ID
  ##        table with columns: "KINASE"      "SUBSTRATE"   "SUB_MOD_RSD" "fdr"         "coef"        "Cancer"      "pair"        "sig"   
  ##        subtype mean of each phosphosite, protein in input table
  cohorts <- subtype_sort
  cancers <- vector(mode = "character", length = length(subtype_sort)); names(cancers) <- subtype_sort
  for (cancer in c("BRCA", "OV", "CO")) {
    tmp <- grepl(cancer, cohorts)
    cancers[tmp] <- cancer
  }
  myvjust = 0.7
  myhjust = 0
  table_bubble$x_print <- paste0(table_bubble$Cancer, "_", "Assoc")
  table_bubble$Cancer_print <- factor(table_bubble$Cancer, levels = c("BRCA", "OV", "CO"))
  y_print <- as.vector(table_bubble$pair)
  mark_pair <- which(table_bubble$shared3can | table_bubble$uniq_BRCA | table_bubble$uniq_OV | table_bubble$uniq_CO)
  y_print[mark_pair] <- paste0("*", y_print[mark_pair])
  table_bubble$y_print <- y_print
  table_bubble_uniq <- unique(table_bubble[,c("KINASE","SUBSTRATE","SUB_MOD_RSD","y_print")])
  
  ## extract phosphosite level
  pho_rsd_split <- formatPhosphosite(table_substrate$Phosphosite, table_substrate$Gene)
  pho_table <- c()
  for (cohort in cohorts) {
    temp <- table_bubble_uniq
    temp$cohort <- cohort
    temp$Cancer <- cancers[cohort]
    temp$pho_subtype <- NA
    for (i in 1:nrow(table_bubble_uniq)) {
      pho_temp <- table_substrate[pho_rsd_split$SUBSTRATE==as.character(temp$SUBSTRATE[i]) & pho_rsd_split$SUB_MOD_RSD==as.character(temp$SUB_MOD_RSD[i]) & !is.na(table_substrate[,cohort]) ,cohort]
      if (length(pho_temp) > 0) {
        temp$pho_subtype[i] <- pho_temp
      }
    }
    pho_table <- rbind(pho_table,temp)
  }
  pho_table <- pho_table[!is.na(pho_table$pho_subtype),]
  
  ## extract kinase protein/phosphoprotein level
  pro_table <- c()
  for (cohort in cohorts) {
    temp <- table_bubble_uniq
    temp$cohort <- cohort
    temp$pro_subtype <- NA
    temp$Cancer <- cancers[cohort]
    for (i in 1:nrow(table_bubble_uniq)) {
      value_tmp <- table_kinase[table_kinase$Gene == as.character(temp$KINASE[i])  & !is.na(table_kinase[,cohort]) ,cohort]
      if (length(value_tmp) > 0){
        temp$pro_subtype[i] <- value_tmp
      }
    }
    pro_table <- rbind(pro_table,temp)
  }
  pro_table <- pro_table[!is.na(pro_table$pro_subtype),]
  
  ## bubble plot
  lim = max(abs(max(table_bubble$coef, na.rm = T)),abs(min(table_bubble$coef, na.rm = T)))
  p = ggplot(table_bubble,aes(x=x_print, y=y_print))# make this the original ethni
  p = p + facet_grid(KINASE~Cancer_print, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p = p + geom_point(aes(color=coef, size =-log10(fdr)),pch=16) 
  p = p + scale_color_gradientn(name= "Coefficient", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
  p = p + theme_bw() + theme_nogrid()
  p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
  p = p + theme(axis.title = element_blank(), axis.text.y = element_text(colour="black", size=10))
  p = p + theme(axis.text.x = element_text(angle = -30, vjust = myvjust, hjust = myhjust))
  p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
  p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
  return(p)
  # plot heatmap for subtype phosphosite level --------------------------
  # lim = max(abs(max(pho_table$pho_subtype, na.rm = T)),abs(min(pho_table$pho_subtype, na.rm = T)))
  # cap <- min(2, lim)
  # pho_table$pho_capped <- as.numeric(pho_table$pho_subtype)
  # pho_table$pho_capped[pho_table$pho_capped > cap] <- cap
  # pho_table$pho_capped[pho_table$pho_capped < (-cap)] <- (-cap)
  # pho_table$cohort_p <- factor(pho_table$cohort, levels = subtype_sort)
  # pho_table$Cancer_print <- factor(pho_table$Cancer, levels = c("BRCA", "OV", "CO"))
  # plot2 = ggplot(pho_table)
  # plot2 = plot2 + geom_tile(aes(x=cohort_p, y=y_print, fill=pho_capped), color=NA)#, linetype="blank") 
  # plot2 = plot2 + facet_grid(KINASE~Cancer_print, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  # plot2 = plot2 + scale_fill_gradientn(name= "phospho_level", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  # plot2 = plot2 + theme_bw() + theme_nogrid()
  # plot2 = plot2 + theme(axis.text.x = element_text(angle = -30, vjust = myvjust, hjust = myhjust))
  # plot2 = plot2 + theme(axis.title = element_blank(),  axis.text.y = element_blank())
  # plot2 = plot2 + theme(axis.ticks=element_blank(),legend.position="bottom")
  # plot2 = plot2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
  # 
  # # plot heatmap for subtype kinase expression level --------------------------
  # lim = max(abs(max(pro_table$pro_subtype, na.rm = T)),abs(min(pro_table$pro_subtype, na.rm = T)))
  # cap <- min(2, lim)
  # pro_table$pro_capped <- pro_table$pro_subtype
  # pro_table$pro_capped[pro_table$pro_capped > cap] <- cap
  # pro_table$pro_capped[pro_table$pro_capped < (-cap)] <- (-cap)
  # pro_table$cohort_p <- factor(pro_table$cohort, levels = subtype_sort)
  # pro_table$Cancer_print <- factor(pro_table$Cancer, levels = c("BRCA", "OV", "CO"))
  # plot3 = ggplot(pro_table)
  # plot3 = plot3 + geom_tile(aes(x=cohort_p, y=y_print, fill=pro_capped), color=NA)#, linetype="blank") 
  # plot3 = plot3 + facet_grid(KINASE~Cancer_print, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  # plot3 = plot3 + scale_fill_gradientn(name= "pro_level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
  # plot3 = plot3 + theme_bw() + theme_nogrid()
  # plot3 = plot3 + theme(axis.title = element_blank(), axis.text.y = element_blank())
  # plot3 = plot3 + theme(axis.text.x = element_text(angle = -30, vjust = myvjust, hjust = myhjust))
  # plot3 = plot3 + theme(axis.ticks=element_blank(), legend.position="bottom")
  # plot3 = plot3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
  # 
  # # plot together --------------------------------------------------------------------
  # fn = paste(outDir, "cptac2p_", cancer_id, "_", protein,'_regrssion_',self, "_", id, '_heatmapWsubtype.pdf',sep ="")
  # grid.newpage()
  # # pdf(fn, height = 14*(len+10)/(top+10), width = 8*(len+10)/(top+10))
  # pdf(fn, height = h, width = w)
  # 
  # plot123 <- cbind(ggplotGrob(p), ggplotGrob(plot2), ggplotGrob(plot3), size = "last")
  # title <- textGrob(paste0(cancer_id, " ", self,"-regulated ",protein, "-substrate pairs(", id ,")"),gp=gpar(fontsize=16))
  # padding <- unit(5,"mm")
  # plottgt <- gtable_add_rows(plot123, 
  #                            heights = grobHeight(title) + padding,
  #                            pos = 0)
  # plottgt <- gtable_add_grob(plottgt, title, 1, 1, 1, ncol(plottgt))
  # grid.draw(plottgt)
  # dev.off()
}

GeneList_bubble_heatmap = function(resultD, inputD, resultDnow, sig_thres, prefix, genelist, id) {
  ks_targets <- vector('list', 2); ks_targets[['cis']] = c('KINASE'); ks_targets[['trans']] <- c('KINASE', 'SUBSTRATE')
  var_list <- c('pro_kin', 'pho_kin')
  names(var_list) <- c('cis', 'trans')
  subtypes_list <- c('pam50', 'MSI_type', 'tumor_normal_WBasal'); names(subtypes_list) <- c("BRCA", "CO", "OV")
  
  for (protein in c('kinase', 'phosphotase')) {
    tn = paste0(resultD,"regression/tables/",protein,"_substrate_regression_cptac2p_3can.txt")
    table_3can <- fread(tn, data.table = F)
    table_3can <- markSigSiteCan(table_3can, sig_thres = sig_thres, protein_type = protein)
    
    for (self in c("cis","trans")) {
      t0 <- table_3can[table_3can$SELF == self & table_3can$fdr_sig & table_3can$coef_sig,]
      var <- var_list[self]
      fdr_var <- paste("FDR_",var,sep = "");  coef_var <- paste("coef_",var,sep = "")
      
      if (nrow(t0) > 0){
        t0 <- t0[order(t0[,fdr_var]), c("KINASE","SUBSTRATE","SUB_MOD_RSD",fdr_var, coef_var,
                                        "Cancer","pair", "shared3can", "uniq_BRCA", "uniq_OV", "uniq_CO")]
        colnames(t0) <- c("KINASE","SUBSTRATE","SUB_MOD_RSD","fdr","coef",
                          "Cancer","pair", "shared3can", "uniq_BRCA", "uniq_OV", "uniq_CO")
        t1 <- t0
        
        ## draw pairs with given kinase/substrate list
        for (ks_target in ks_targets[[self]]) {
          table_bubble <- t1[t1[,ks_target] %in% genelist, ]
          if (nrow(table_bubble) > 0){
            table_substrate <- NULL
            table_kinase <- NULL
            for (cancer in c("CO", "BRCA", "OV")) {
              pho_subtype_mean <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_",subtypes_list[cancer],"_mean.txt"), data.table = F)
              pro_subtype_mean <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_",subtypes_list[cancer],"_mean.txt"), data.table = F)
              phog_subtype_mean <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_collapsed_",subtypes_list[cancer],"_mean.txt"), data.table = F)
              if (cancer == "OV") {
                pho_subtype_mean <- pho_subtype_mean[, c("Gene", "Phosphosite","Tumor", "Normal")]
                pro_subtype_mean <- pro_subtype_mean[, c("Gene", "Tumor", "Normal")]
                phog_subtype_mean <- phog_subtype_mean[, c("Gene","Tumor", "Normal")]
              }
              if (self == "cis") {
                table_kinase_tmp <- pro_subtype_mean
              } else {
                table_kinase_tmp <- phog_subtype_mean
              }
              table_substrate_tmp <- pho_subtype_mean
              ## edit subtype names to include cancer types
              colname_tmp <- colnames(table_kinase_tmp)[!(colnames(table_kinase_tmp) %in% c("Gene", "Phosphosite"))]
              colnames(table_kinase_tmp)[!(colnames(table_kinase_tmp) %in% c("Gene", "Phosphosite"))] <- paste0(cancer, "_", colname_tmp)
              colname_tmp <- colnames(table_substrate_tmp)[!(colnames(table_substrate_tmp) %in% c("Gene", "Phosphosite"))]
              colnames(table_substrate_tmp)[!(colnames(table_substrate_tmp) %in% c("Gene", "Phosphosite"))] <- paste0(cancer, "_", colname_tmp)
              table_kinase_tmp$Cancer <- cancer
              table_substrate_tmp$Cancer <- cancer
              if (is.null(table_substrate)) {
                table_substrate <- table_substrate_tmp[table_substrate_tmp$Gene %in% table_bubble$SUBSTRATE,]
                table_kinase <- table_kinase_tmp[table_kinase_tmp$Gene %in% table_bubble$KINASE,]
              } else {
                table_substrate <- merge(table_substrate, table_substrate_tmp, by = c("Gene", "Phosphosite", "Cancer"), all = T)
                table_kinase <- merge(table_kinase, table_kinase_tmp, by = c("Gene", "Cancer"), all = T)
              }
            }
            bubble_heatmap("3can", protein, self, 
                           paste0(id, ks_target), table_bubble, table_substrate, table_kinase, resultDnow, 
                           subtype_sort = c("BRCA_LumA", "BRCA_LumB", "BRCA_Her2", "BRCA_Basal" ,"BRCA_Normal", "OV_Tumor","OV_Normal","CO_MSI","CO_MSS","CO_other","CO_Normal"),
                           h = 15, w = 8)
          } else {
            print('not enough data')
          }
        }
      }
    }
  }
}

kspairList_bubble_heatmap = function(resultD, inputD, resultDnow, sig_thres, prefix, pairlist, id, table_3can) {
  var_list <- c('pro_kin', 'pho_kin')
  names(var_list) <- c('cis', 'trans')
  subtypes_list <- c('pam50', 'MSI_type', 'tumor_normal_WBasal'); names(subtypes_list) <- c("BRCA", "CO", "OV")
  
  for (protein in unique(table_3can$enzyme_type)) {
    # tn = paste0(resultD,"regression/tables/",protein,"_substrate_regression_cptac2p_3can.txt")
    # table_3can <- fread(tn, data.table = F)
    if (!("uniq_BRCA" %in% colnames(table_3can))) {
      table_3can <- markSigSiteCan(table_3can, sig_thres = sig_thres, protein_type = protein)
    }
    
    for (self in c("cis","trans")) {
      t0 <- table_3can[table_3can$SELF == self & table_3can$fdr_sig & table_3can$coef_sig,]
      var <- var_list[self]
      fdr_var <- paste("FDR_",var,sep = "");  coef_var <- paste("coef_",var,sep = "")
      
      if (nrow(t0) > 0){
        t0 <- t0[order(t0[,fdr_var]), c("KINASE","SUBSTRATE","SUB_MOD_RSD",fdr_var, coef_var,
                                        "Cancer","pair", "shared3can", "uniq_BRCA", "uniq_OV", "uniq_CO")]
        colnames(t0) <- c("KINASE","SUBSTRATE","SUB_MOD_RSD","fdr","coef",
                          "Cancer","pair", "shared3can", "uniq_BRCA", "uniq_OV", "uniq_CO")
        t1 <- t0
        
        ## draw pairs with given kinase/substrate list
        table_bubble <- t1[t1$pair %in% pairlist, ]
        if (nrow(table_bubble) > 0){
          table_substrate <- NULL
          table_kinase <- NULL
          for (cancer in c("CO", "BRCA", "OV")) {
            pho_subtype_mean <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_",subtypes_list[cancer],"_mean.txt"), data.table = F)
            pro_subtype_mean <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_",subtypes_list[cancer],"_mean.txt"), data.table = F)
            phog_subtype_mean <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_collapsed_",subtypes_list[cancer],"_mean.txt"), data.table = F)
            if (cancer == "OV") {
              pho_subtype_mean <- pho_subtype_mean[, c("Gene", "Phosphosite","Tumor", "Normal")]
              pro_subtype_mean <- pro_subtype_mean[, c("Gene", "Tumor", "Normal")]
              phog_subtype_mean <- phog_subtype_mean[, c("Gene","Tumor", "Normal")]
            }
            if (self == "cis") {
              table_kinase_tmp <- pro_subtype_mean
            } else {
              table_kinase_tmp <- phog_subtype_mean
            }
            table_substrate_tmp <- pho_subtype_mean
            ## edit subtype names to include cancer types
            colname_tmp <- colnames(table_kinase_tmp)[!(colnames(table_kinase_tmp) %in% c("Gene", "Phosphosite"))]
            colnames(table_kinase_tmp)[!(colnames(table_kinase_tmp) %in% c("Gene", "Phosphosite"))] <- paste0(cancer, "_", colname_tmp)
            colname_tmp <- colnames(table_substrate_tmp)[!(colnames(table_substrate_tmp) %in% c("Gene", "Phosphosite"))]
            colnames(table_substrate_tmp)[!(colnames(table_substrate_tmp) %in% c("Gene", "Phosphosite"))] <- paste0(cancer, "_", colname_tmp)
            table_kinase_tmp$Cancer <- cancer
            table_substrate_tmp$Cancer <- cancer
            if (is.null(table_substrate)) {
              table_substrate <- table_substrate_tmp[table_substrate_tmp$Gene %in% table_bubble$SUBSTRATE,]
              table_kinase <- table_kinase_tmp[table_kinase_tmp$Gene %in% table_bubble$KINASE,]
            } else {
              table_substrate <- merge(table_substrate, table_substrate_tmp, by = c("Gene", "Phosphosite", "Cancer"), all = T)
              table_kinase <- merge(table_kinase, table_kinase_tmp, by = c("Gene", "Cancer"), all = T)
            }
          }
          p <- bubble_heatmap("3can", protein, self, 
                         paste0(id), table_bubble, table_substrate, table_kinase, resultDnow, 
                         subtype_sort = c("BRCA_LumA", "BRCA_LumB", "BRCA_Her2", "BRCA_Basal" ,"BRCA_Normal", "OV_Tumor","OV_Normal","CO_MSI","CO_MSS","CO_other","CO_Normal"),
                         h = 20, w = 8)
          return(p)
        } else {
          print('not enough data')
        }
      }
    }
  }
}
