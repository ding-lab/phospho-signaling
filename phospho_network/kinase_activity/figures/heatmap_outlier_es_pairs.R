# Yige Wu @ WashU 2018 Mar
# annotate table f hot and cold kinase-substrate pairs by fisher's test

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# set variables -----------------------------------------------------------
var_list <- c('pro_kin', 'pho_kin')
names(var_list) <- c('cis', 'trans')
subtypes_list <- c('pam50', 'MSI_type', 'tumor_normal_WBasal'); names(subtypes_list) <- c("BRCA", "CO", "OV")
subtype_sort = c("BRCA_LumA", "BRCA_LumB", "BRCA_Her2", "BRCA_Basal" ,"BRCA_Normal", "OV_Tumor","OV_Normal","CO_MSI","CO_MSS","CO_other","CO_Normal")
myvjust = 0.7
myhjust = 0
cancers_sort <- c("BRCA", "OV", "CO")

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
sup_cans_tab_pk <- markSigSiteCan(sup_cans_tab_pk, sig_thres = 0.05, enzyme_type = enzyme_type)


# drivers for all cancers -------------------------------------------------
driver_cans <- NULL
for (cancer in cancers_sort) {
  tmp <- loadGeneList(gene_type = "driver", cancer = cancer, is.soft.limit = "soft")
  driver_cans <- c(driver_cans, tmp)
}

driver_cans <- unique(driver_cans)
pairlist <- NULL
for (cancer in cancers_sort) {
  tmp <- unique(sup_cans_tab_pk$pair[!is.na(sup_cans_tab_pk[, paste0("uniq_", cancer)]) & sup_cans_tab_pk[, paste0("uniq_", cancer)]  & sup_cans_tab_pk$Cancer == cancer & (sup_cans_tab_pk$GENE %in% driver_cans | sup_cans_tab_pk$SUB_GENE %in% driver_cans)])
  pairlist <- unique(c(pairlist, tmp))
}

tab_cans_tmp <- sup_cans_tab_pk[sup_cans_tab_pk$pair %in% pairlist,]

for (self in c("trans")) {
  t0 <- tab_cans_tmp[tab_cans_tmp$SELF == self,]
  var <- var_list[self]
  fdr_var <- paste("FDR_",var,sep = "");  coef_var <- paste("coef_",var,sep = "")
  
  if (nrow(t0) > 0){
    t0 <- t0[order(t0[,fdr_var]), c("GENE","SUB_GENE","SUB_MOD_RSD",fdr_var, coef_var,
                                    "Cancer","pair", "shared3can", "uniq_BRCA", "uniq_OV", "uniq_CO")]
    colnames(t0) <- c("KINASE","SUBSTRATE","SUB_MOD_RSD","fdr","coef",
                      "Cancer","pair", "shared3can", "uniq_BRCA", "uniq_OV", "uniq_CO")
    t1 <- t0
    
    ## draw pairs with given kinase/substrate list
    table_bubble <- t1[t1$pair %in% pairlist, ]
    if (nrow(table_bubble) > 0){
      table_substrate <- NULL
      table_kinase <- NULL
      for (cancer2 in unique(tab_cans_tmp$Cancer)) {
        pho_subtype_mean <- fread(input = paste0(cptac_sharedD, cancer2, "/", prefix[cancer2], "_PHO_formatted_normalized_",subtypes_list[cancer2],"_mean.txt"), data.table = F)
        pro_subtype_mean <- fread(input = paste0(cptac_sharedD, cancer2, "/", prefix[cancer2], "_PRO_formatted_normalized_",subtypes_list[cancer2],"_mean.txt"), data.table = F)
        phog_subtype_mean <- fread(input = paste0(cptac_sharedD, cancer2, "/", prefix[cancer2], "_PHO_formatted_normalized_collapsed_",subtypes_list[cancer2],"_mean.txt"), data.table = F)
        if (cancer2 == "OV") {
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
        colnames(table_kinase_tmp)[!(colnames(table_kinase_tmp) %in% c("Gene", "Phosphosite"))] <- paste0(cancer2, "_", colname_tmp)
        colname_tmp <- colnames(table_substrate_tmp)[!(colnames(table_substrate_tmp) %in% c("Gene", "Phosphosite"))]
        colnames(table_substrate_tmp)[!(colnames(table_substrate_tmp) %in% c("Gene", "Phosphosite"))] <- paste0(cancer2, "_", colname_tmp)
        table_kinase_tmp$Cancer <- cancer2
        table_substrate_tmp$Cancer <- cancer2
        if (is.null(table_substrate)) {
          table_substrate <- table_substrate_tmp[table_substrate_tmp$Gene %in% table_bubble$SUBSTRATE,]
          table_kinase <- table_kinase_tmp[table_kinase_tmp$Gene %in% table_bubble$KINASE,]
        } else {
          table_substrate <- merge(table_substrate, table_substrate_tmp, by = c("Gene", "Phosphosite", "Cancer"), all = T)
          table_kinase <- merge(table_kinase, table_kinase_tmp, by = c("Gene", "Cancer"), all = T)
        }
      }
      
      cohorts <- subtype_sort
      cancers <- vector(mode = "character", length = length(subtype_sort)); names(cancers) <- subtype_sort
      for (cancer2 in unique(tab_cans_tmp$Cancer)) {
        tmp <- grepl(cancer2, cohorts)
        cancers[tmp] <- cancer2
      }
      
      table_bubble$x_print <- paste0(table_bubble$Cancer, "_", "Assoc")
      table_bubble$Cancer_print <- factor(table_bubble$Cancer, levels = cancers_sort)
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
      # lim = max(abs(max(table_bubble$coef, na.rm = T)),abs(min(table_bubble$coef, na.rm = T)))
      coef_lim <- 1
      table_bubble$coef_capped <- table_bubble$coef
      table_bubble$coef_capped[table_bubble$coef_capped > coef_lim] <- coef_lim
      table_bubble$coef_capped[table_bubble$coef_capped < (-coef_lim)] <- (-coef_lim)
      
      ## mark cancer unique ones
      star <- vector(mode = "logical", length = nrow(table_bubble))
      for (cancer2 in cancers_sort) {
        star[table_bubble$Cancer == cancer2 & table_bubble[, paste0('uniq_', cancer2)]] <- TRUE
      }
      table_bubble$star <- star
      p = ggplot(table_bubble,aes(x=x_print, y=y_print))# make this the original ethni
      p = p + geom_point(aes(color=coef_capped, size =-log10(fdr)),pch=16) 
      p = p + geom_point(data = table_bubble[table_bubble$star,], mapping = aes(x=x_print, y=y_print),  shape = 8, size = 2) 
      p = p + scale_color_gradientn(name= "Coefficient", na.value=NA, colours=RdBu1024, limit=c(-coef_lim,coef_lim))
      p = p + facet_grid(KINASE~Cancer_print, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
      p = p + theme_bw() + theme_nogrid()
      p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
      p = p + theme(axis.title = element_blank(), axis.text.y = element_text(colour="black", size=10))
      p = p + theme(axis.text.x = element_text(angle = -30, vjust = myvjust, hjust = myhjust))
      p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
      p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
      
      # plot heatmap for subtype phosphosite level --------------------------
      lim = max(abs(max(pho_table$pho_subtype, na.rm = T)),abs(min(pho_table$pho_subtype, na.rm = T)))
      cap <- min(2, lim)
      pho_table$pho_capped <- as.numeric(pho_table$pho_subtype)
      pho_table$pho_capped[pho_table$pho_capped > cap] <- cap
      pho_table$pho_capped[pho_table$pho_capped < (-cap)] <- (-cap)
      pho_table$cohort_p <- factor(pho_table$cohort, levels = subtype_sort)
      pho_table$Cancer_print <- factor(pho_table$Cancer, levels = c("BRCA", "OV", "CO"))
      plot2 = ggplot(pho_table)
      plot2 = plot2 + geom_tile(aes(x=cohort_p, y=y_print, fill=pho_capped), color=NA)#, linetype="blank")
      plot2 = plot2 + facet_grid(KINASE~Cancer_print, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
      plot2 = plot2 + scale_fill_gradientn(name= "phospho_level", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
      plot2 = plot2 + theme_bw() + theme_nogrid()
      plot2 = plot2 + theme(axis.text.x = element_text(angle = -30, vjust = myvjust, hjust = myhjust))
      plot2 = plot2 + theme(axis.title = element_blank(),  axis.text.y = element_blank())
      plot2 = plot2 + theme(axis.ticks=element_blank(),legend.position="bottom")
      plot2 = plot2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
      
      # plot heatmap for subtype kinase expression level --------------------------
      lim = max(abs(max(pro_table$pro_subtype, na.rm = T)),abs(min(pro_table$pro_subtype, na.rm = T)))
      cap <- min(2, lim)
      pro_table$pro_capped <- pro_table$pro_subtype
      pro_table$pro_capped[pro_table$pro_capped > cap] <- cap
      pro_table$pro_capped[pro_table$pro_capped < (-cap)] <- (-cap)
      pro_table$cohort_p <- factor(pro_table$cohort, levels = subtype_sort)
      pro_table$Cancer_print <- factor(pro_table$Cancer, levels = c("BRCA", "OV", "CO"))
      plot3 = ggplot(pro_table)
      plot3 = plot3 + geom_tile(aes(x=cohort_p, y=y_print, fill=pro_capped), color=NA)#, linetype="blank")
      plot3 = plot3 + facet_grid(KINASE~Cancer_print, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
      plot3 = plot3 + scale_fill_gradientn(name= "pro_level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
      plot3 = plot3 + theme_bw() + theme_nogrid()
      plot3 = plot3 + theme(axis.title = element_blank(), axis.text.y = element_blank())
      plot3 = plot3 + theme(axis.text.x = element_text(angle = -30, vjust = myvjust, hjust = myhjust))
      plot3 = plot3 + theme(axis.ticks=element_blank(), legend.position="bottom")
      plot3 = plot3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
      
      # plot together --------------------------------------------------------------------
      id <- paste0('uniq_cancer')
      fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_regrssion_',self, "_", paste0(cancers_sort, collapse = '&') , 'uniq_heatmapWsubtype.pdf',sep ="")
      grid.newpage()
      # pdf(fn, height = 14*(len+10)/(top+10), width = 8*(len+10)/(top+10))
      pdf(fn, height = 20, width = 8, useDingbats = FALSE)
      
      plot123 <- cbind(ggplotGrob(p), ggplotGrob(plot2), ggplotGrob(plot3), size = "last")
      title <- textGrob(paste0(paste0(cancers_sort, collapse = '&'), " ", self,"-regulated ", enzyme_type, "-substrate pairs(", id ,")"),gp=gpar(fontsize=16))
      padding <- unit(5,"mm")
      plottgt <- gtable_add_rows(plot123,
                                 heights = grobHeight(title) + padding,
                                 pos = 0)
      plottgt <- gtable_add_grob(plottgt, title, 1, 1, 1, ncol(plottgt))
      grid.draw(plottgt)
      dev.off()
    } else {
      print('not enough data')
    }
  }
}



# plot each set cancer drivers --------------------------------------------
for (cancer in cancers_sort) {
  # pairlist <- unique(sup_cans_tab_pk$pair[sup_cans_tab_pk[, paste0("uniq_", cancer)]  & sup_cans_tab_pk$Cancer == cancer & sup_cans_tab_pk$fisher_pvalues < 0.1 & sup_cans_tab_pk$fisher_pvalues > 0])
  pairlist <- unique(sup_cans_tab_pk$pair[!is.na(sup_cans_tab_pk[, paste0("uniq_", cancer)]) & sup_cans_tab_pk[, paste0("uniq_", cancer)]  & sup_cans_tab_pk$Cancer == cancer & (sup_cans_tab_pk$enzyme.driver_type != "notdriver" | sup_cans_tab_pk$substrate.driver_type != "notdriver") & !is.na(sup_cans_tab_pk$substrate.driver_type)])
  # pairlist <- unique(sup_cans_tab_pk$pair[!is.na(sup_cans_tab_pk[, paste0("uniq_", cancer)]) & sup_cans_tab_pk[, paste0("uniq_", cancer)]  & sup_cans_tab_pk$Cancer == cancer & !is.na(sup_cans_tab_pk$substrate.driver_type)])
  # pairlist <- unique(sup_cans_tab_pk$pair[!is.na(sup_cans_tab_pk[, paste0("uniq_", cancer)]) & sup_cans_tab_pk$Cancer == cancer & (sup_cans_tab_pk$enzyme.driver_type != "notdriver" | sup_cans_tab_pk$substrate.driver_type != "notdriver") & !is.na(sup_cans_tab_pk$substrate.driver_type)])
  # 
  tab_cans_tmp <- sup_cans_tab_pk[sup_cans_tab_pk$pair %in% pairlist,]
  # kins_tab <- data.frame(table(tab_cans_tmp$GENE))
  # kins_tab <- kins_tab[order(kins_tab$Freq, decreasing = T),]
  for (self in c("cis","trans")) {
    t0 <- tab_cans_tmp[tab_cans_tmp$SELF == self,]
    var <- var_list[self]
    fdr_var <- paste("FDR_",var,sep = "");  coef_var <- paste("coef_",var,sep = "")
    
    if (nrow(t0) > 0){
      t0 <- t0[order(t0[,fdr_var]), c("GENE","SUB_GENE","SUB_MOD_RSD",fdr_var, coef_var,
                                      "Cancer","pair", "shared3can", "uniq_BRCA", "uniq_OV", "uniq_CO")]
      colnames(t0) <- c("KINASE","SUBSTRATE","SUB_MOD_RSD","fdr","coef",
                        "Cancer","pair", "shared3can", "uniq_BRCA", "uniq_OV", "uniq_CO")
      t1 <- t0
      
      ## draw pairs with given kinase/substrate list
      table_bubble <- t1[t1$pair %in% pairlist, ]
      if (nrow(table_bubble) > 0){
        table_substrate <- NULL
        table_kinase <- NULL
        for (cancer2 in unique(tab_cans_tmp$Cancer)) {
          pho_subtype_mean <- fread(input = paste0(cptac_sharedD, cancer2, "/", prefix[cancer2], "_PHO_formatted_normalized_",subtypes_list[cancer2],"_mean.txt"), data.table = F)
          pro_subtype_mean <- fread(input = paste0(cptac_sharedD, cancer2, "/", prefix[cancer2], "_PRO_formatted_normalized_",subtypes_list[cancer2],"_mean.txt"), data.table = F)
          phog_subtype_mean <- fread(input = paste0(cptac_sharedD, cancer2, "/", prefix[cancer2], "_PHO_formatted_normalized_collapsed_",subtypes_list[cancer2],"_mean.txt"), data.table = F)
          if (cancer2 == "OV") {
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
          colnames(table_kinase_tmp)[!(colnames(table_kinase_tmp) %in% c("Gene", "Phosphosite"))] <- paste0(cancer2, "_", colname_tmp)
          colname_tmp <- colnames(table_substrate_tmp)[!(colnames(table_substrate_tmp) %in% c("Gene", "Phosphosite"))]
          colnames(table_substrate_tmp)[!(colnames(table_substrate_tmp) %in% c("Gene", "Phosphosite"))] <- paste0(cancer2, "_", colname_tmp)
          table_kinase_tmp$Cancer <- cancer2
          table_substrate_tmp$Cancer <- cancer2
          if (is.null(table_substrate)) {
            table_substrate <- table_substrate_tmp[table_substrate_tmp$Gene %in% table_bubble$SUBSTRATE,]
            table_kinase <- table_kinase_tmp[table_kinase_tmp$Gene %in% table_bubble$KINASE,]
          } else {
            table_substrate <- merge(table_substrate, table_substrate_tmp, by = c("Gene", "Phosphosite", "Cancer"), all = T)
            table_kinase <- merge(table_kinase, table_kinase_tmp, by = c("Gene", "Cancer"), all = T)
          }
        }
        
        cohorts <- subtype_sort
        cancers <- vector(mode = "character", length = length(subtype_sort)); names(cancers) <- subtype_sort
        for (cancer2 in unique(tab_cans_tmp$Cancer)) {
          tmp <- grepl(cancer2, cohorts)
          cancers[tmp] <- cancer2
        }

        table_bubble$x_print <- paste0(table_bubble$Cancer, "_", "Assoc")
        table_bubble$Cancer_print <- factor(table_bubble$Cancer, levels = cancers_sort)
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
        # lim = max(abs(max(table_bubble$coef, na.rm = T)),abs(min(table_bubble$coef, na.rm = T)))
        coef_lim <- 1
        table_bubble$coef_capped <- table_bubble$coef
        table_bubble$coef_capped[table_bubble$coef_capped > coef_lim] <- coef_lim
        table_bubble$coef_capped[table_bubble$coef_capped < (-coef_lim)] <- (-coef_lim)
        
        p = ggplot(table_bubble,aes(x=x_print, y=y_print))# make this the original ethni
        p = p + facet_grid(KINASE~Cancer_print, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
        p = p + geom_point(aes(color=coef_capped, size =-log10(fdr)),pch=16) 
        p = p + scale_color_gradientn(name= "Coefficient", na.value=NA, colours=RdBu1024, limit=c(-coef_lim,coef_lim))
        p = p + theme_bw() + theme_nogrid()
        p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
        p = p + theme(axis.title = element_blank(), axis.text.y = element_text(colour="black", size=10))
        p = p + theme(axis.text.x = element_text(angle = -30, vjust = myvjust, hjust = myhjust))
        p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
        p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
        
        # plot heatmap for subtype phosphosite level --------------------------
        lim = max(abs(max(pho_table$pho_subtype, na.rm = T)),abs(min(pho_table$pho_subtype, na.rm = T)))
        cap <- min(2, lim)
        pho_table$pho_capped <- as.numeric(pho_table$pho_subtype)
        pho_table$pho_capped[pho_table$pho_capped > cap] <- cap
        pho_table$pho_capped[pho_table$pho_capped < (-cap)] <- (-cap)
        pho_table$cohort_p <- factor(pho_table$cohort, levels = subtype_sort)
        pho_table$Cancer_print <- factor(pho_table$Cancer, levels = c("BRCA", "OV", "CO"))
        plot2 = ggplot(pho_table)
        plot2 = plot2 + geom_tile(aes(x=cohort_p, y=y_print, fill=pho_capped), color=NA)#, linetype="blank")
        plot2 = plot2 + facet_grid(KINASE~Cancer_print, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
        plot2 = plot2 + scale_fill_gradientn(name= "phospho_level", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
        plot2 = plot2 + theme_bw() + theme_nogrid()
        plot2 = plot2 + theme(axis.text.x = element_text(angle = -30, vjust = myvjust, hjust = myhjust))
        plot2 = plot2 + theme(axis.title = element_blank(),  axis.text.y = element_blank())
        plot2 = plot2 + theme(axis.ticks=element_blank(),legend.position="bottom")
        plot2 = plot2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))

        # plot heatmap for subtype kinase expression level --------------------------
        lim = max(abs(max(pro_table$pro_subtype, na.rm = T)),abs(min(pro_table$pro_subtype, na.rm = T)))
        cap <- min(2, lim)
        pro_table$pro_capped <- pro_table$pro_subtype
        pro_table$pro_capped[pro_table$pro_capped > cap] <- cap
        pro_table$pro_capped[pro_table$pro_capped < (-cap)] <- (-cap)
        pro_table$cohort_p <- factor(pro_table$cohort, levels = subtype_sort)
        pro_table$Cancer_print <- factor(pro_table$Cancer, levels = c("BRCA", "OV", "CO"))
        plot3 = ggplot(pro_table)
        plot3 = plot3 + geom_tile(aes(x=cohort_p, y=y_print, fill=pro_capped), color=NA)#, linetype="blank")
        plot3 = plot3 + facet_grid(KINASE~Cancer_print, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
        plot3 = plot3 + scale_fill_gradientn(name= "pro_level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
        plot3 = plot3 + theme_bw() + theme_nogrid()
        plot3 = plot3 + theme(axis.title = element_blank(), axis.text.y = element_blank())
        plot3 = plot3 + theme(axis.text.x = element_text(angle = -30, vjust = myvjust, hjust = myhjust))
        plot3 = plot3 + theme(axis.ticks=element_blank(), legend.position="bottom")
        plot3 = plot3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))

        # plot together --------------------------------------------------------------------
        id <- paste0('uniq_', cancer)
        fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_regrssion_',self, "_", cancer , 'uniq_heatmapWsubtype.pdf',sep ="")
        grid.newpage()
        # pdf(fn, height = 14*(len+10)/(top+10), width = 8*(len+10)/(top+10))
        pdf(fn, height = 15, width = 8)

        plot123 <- cbind(ggplotGrob(p), ggplotGrob(plot2), ggplotGrob(plot3), size = "last")
        title <- textGrob(paste0(cancer, " ", self,"-regulated ", enzyme_type, "-substrate pairs(", id ,")"),gp=gpar(fontsize=16))
        padding <- unit(5,"mm")
        plottgt <- gtable_add_rows(plot123,
                                   heights = grobHeight(title) + padding,
                                   pos = 0)
        plottgt <- gtable_add_grob(plottgt, title, 1, 1, 1, ncol(plottgt))
        grid.draw(plottgt)
        dev.off()
      } else {
        print('not enough data')
      }
    }
  }
}
