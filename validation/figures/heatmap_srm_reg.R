# Yige Wu @ WashU 2018 Apr
## make heatmap showing the ks outlier samples really agree regarding the order of abundance


# source ------------------------------------------------------------------
library(readxl)
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# inputs ------------------------------------------------------------------
## input the date of receiving the data
email_date <- "2018-04-20"

## input the sample ID mapping info
clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180307.txt"), sep = "\t")

## input the SRM data
# cis_tab <- read_excel("~/Box Sync/cptac2p/tao/2018-03-29/Druggable-kinase-outlier-CPTAC2-SRM_03292018_omnipath_result_commented.xlsx",
#                       sheet = "cis")
cis_tab <- read_excel("~/Box Sync/cptac2p/tao/2018-04-20/Druggable-kinase-outlier-CPTAC2prospective-SRM_04192018.xlsx",
                      sheet = "cis")
cis_tab <- data.frame(cis_tab)
# trans_tab <- read_excel("~/Box Sync/cptac2p/tao/2018-03-29/Druggable-kinase-outlier-CPTAC2-SRM_03292018_omnipath_result_commented.xlsx",
#                       sheet = "trans")
trans_tab <- read_excel("~/Box Sync/cptac2p/tao/2018-04-20/Druggable-kinase-outlier-CPTAC2prospective-SRM_04192018.xlsx",
                      sheet = "trans")
trans_tab <- data.frame(trans_tab)
pp_tab <- read_excel(path = paste0("~/Box Sync/cptac2p/tao/", email_date, "/Phosphopeptide_Crosstab.xlsx"))
pp_tab <- data.frame(pp_tab)
## annotate the cis and trans tab to examine if there's mismatch of particiapnt IDs
test <- merge(cis_tab, clinical[,1:5], by.x = c("Specimen_Label"), by.y = c("Specimen.Label"), all.x = T, sort = F)
test <- merge(trans_tab, clinical[,1:5], by.x = c("Specimen_Label"), by.y = c("Specimen.Label"), all.x = T, sort = F)



# fill in the blanks ------------------------------------------------------
tab <- cis_tab
for (cols in c(colnames(tab)[1:7])) {
  tmp <- as.vector(tab[, cols])
  for (i in which(is.na(tmp))) {
    tmp[i] <- tmp[i-1]
  }
  tab[, cols] <- tmp
}
cis_tab <- tab

tab <- trans_tab
for (cols in c(colnames(tab)[1:7])) {
  tmp <- as.vector(tab[, cols])
  for (i in which(is.na(tmp))) {
    tmp[i] <- tmp[i-1]
  }
  tab[, cols] <- tmp
}
trans_tab <- tab

# CO all samples --------------------------------------------------------------------
colnames(trans_tab)[5] <- "peptide"
common_columns <- intersect(colnames(trans_tab), colnames(cis_tab))
tab <- rbind(trans_tab[, common_columns], cis_tab[, common_columns])
for (cancer in c("CO")) {
  pp_cancer <- pp_tab[, c("Peptide.Sequence", "Protein.Name", colnames(pp_tab)[grepl(pattern = cancer, x = colnames(pp_tab))])]
  tab_cancer <- tab[tab$cancer_type ==  cancer,]
  samples <- unique(tab_cancer[, c("Specimen_Label", "Participant_ID")])
  samples <- merge(samples, clinical[, c("Specimen.Label", "Participant.ID", "tumor_normal")], by.x = c("Specimen_Label"), by.y = c("Specimen.Label"), all.x = T)
  samples_dup <- unique(samples$Participant_ID[duplicated(samples$Participant_ID)])
  samples_nondup <- samples[!(samples %in% samples_dup),]
  
  ## input global data
  phog <- fread(input = paste(inputD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged.txt",sep=""), data.table = F)
  
  pho <- fread(input = paste(inputD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_replicate_averaged.txt",sep=""), data.table = F)
  pho <- cbind(formatPhosphosite(phosphosite_vector = pho$Phosphosite, pho$Gene), pho[,!(colnames(pho) %in% c("Gene", "Phosphosite"))])
  
  for (peptide in unique(pp_cancer$Peptide.Sequence)) {
    pp_cancer_pep <- pp_cancer[pp_cancer$Peptide.Sequence == toupper(peptide),]
    tab_cancer_pep <- tab_cancer[grepl(pattern = peptide, x = tab_cancer$peptide, ignore.case = T),]
    pp_cancer_pep.m <- melt(pp_cancer_pep, id.vars = c("Peptide.Sequence", "Protein.Name"))
    pp_cancer_pep.m$Participant_ID <- str_split_fixed(string = pp_cancer_pep.m$variable, pattern = "X", 2)[,2]
    pp_cancer_pep.m.nodup <- pp_cancer_pep.m[!(pp_cancer_pep.m$Participant_ID %in% samples_dup),]
    pp_cancer_pep.m.nodup <- merge(pp_cancer_pep.m.nodup, samples, by = c("Participant_ID"), all.x = T)
    
    sub <- unique(tab_cancer_pep$SUBSTRATE)
    rsd <- unique(tab_cancer_pep$SUB_MOD_RSD)
    pho_sub <- pho[pho$SUBSTRATE == sub & pho$SUB_MOD_RSD == rsd,]
    if (nrow(pho_sub) > 0 ){
      pho_sub.m <- melt(pho_sub, id.vars = c("SUBSTRATE", "transcript", "SUB_MOD_RSD"))
      colnames(pho_sub.m) <- c("substrate", "transcript", "SUB_MOD_RSD", "sample", "pho_sub")
      kinases <- unique(tab_cancer_pep$KINASE)
      for (kinase in kinases) {
        phog_kin <- phog[phog$Gene == kinase,]
        phog_kin.m <- melt(phog_kin, id.vars = c("Gene"))
        colnames(phog_kin.m) <- c("kinase", "sample", "phog_kin")
        df <- merge(pho_sub.m, phog_kin.m, all = T)
        df <- merge(df, pp_cancer_pep.m.nodup, by.x = c("sample"), by.y = c("Specimen_Label"), all = T)
        
        df$tumor_normal <- NULL
        df <- merge(df, clinical[, c("Specimen.Label", "tumor_normal")], by.x = c("sample"), by.y = c("Specimen.Label"), all.x = T)
        df <- df[!is.na(df$tumor_normal),]
        df$kinase <- kinase
        df <- merge(df, tab_cancer_pep, by.x = c("sample"), by.y = c("Specimen_Label"), all.x = T)
        df <- df[!is.na(df$phog_kin)||!is.na(df$pho_sub),]
        df$sample <- as.vector(df$sample)
        df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$phog_kin)])
        df$tumor_normal <- factor(df$tumor_normal, levels = c("Tumor", "Normal"))
        
        p0 = ggplot(df)
        p0 = p0 + geom_tile(aes(x=sample, y= "Tumor (purple); Normal (green)", fill= tumor_normal, color = outlier_type, width=0.7, height=0.7), size=0.5)
        p0 = p0 + scale_fill_manual(values = c("Tumor" = "#6A3D9A", "Normal" = "#33A02C"))
        p0 = p0 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
        p0 = p0 + theme_bw() + theme_nogrid()
        p0 = p0 + theme(axis.text.x = element_blank())
        p0 = p0 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                        strip.background = element_blank())
        p0 = p0 + theme(axis.ticks=element_blank(), 
                        legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
        p0 = p0 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
        p0 = p0 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
        
        cap <- min(max(abs(df$phog_kin), na.rm = T),2)
        df$phog_kin_capped <- df$phog_kin
        df$phog_kin_capped[df$phog_kin > cap] <- cap
        df$phog_kin_capped[df$phog_kin < (-cap)] <- (-cap)
        p1 = ggplot(df)
        p1 = p1 + geom_tile(aes(x=sample, y= paste0(kinase, " phospho (global)"), fill= phog_kin_capped, color = outlier_type, width=0.7, height=0.7), size=0.5)
        p1 = p1 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
        p1 = p1 + scale_fill_gradientn(name= "global kinase phosphorylation abundance", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
        p1 = p1 + theme_bw() + theme_nogrid()
        p1 = p1 + theme(axis.text.x = element_blank())
        p1 = p1 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                        strip.background = element_blank())
        p1 = p1 + theme(axis.ticks=element_blank(), 
                        legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
        p1 = p1 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
        p1 = p1 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
        
        cap <- min(max(abs(df$pho_sub), na.rm = T),2)
        df$pho_sub_capped <- df$pho_sub
        df$pho_sub_capped[df$pho_sub > cap] <- cap
        df$pho_sub_capped[df$pho_sub < (-cap)] <- (-cap)
        p2 = ggplot(df)
        p2 = p2 + geom_tile(aes(x=sample, y= paste0(sub, ":", rsd, "(global)"), fill= pho_sub_capped, color = outlier_type, width=0.7, height=0.7), size=0.5)
        p2 = p2 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
        p2 = p2 + scale_fill_gradientn(name= "global substrate phosphosite abundance", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
        p2 = p2 + theme_bw() + theme_nogrid()
        p2 = p2 + theme(axis.text.x = element_blank())
        p2 = p2 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                        strip.background = element_blank(), strip.text.x = element_blank())
        p2 = p2 + theme(axis.ticks=element_blank(), 
                        legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
        p2 = p2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
        p2 = p2 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
        
        df$value[df$value == "#N/A"] <- NA
        df$value <- as.numeric(as.vector(df$value))
        df$value[is.na(df$value) & !is.na(df$phosphopeptide.heavy.ratio)] <- as.numeric(df$phosphopeptide.heavy.ratio[is.na(df$value) & !is.na(df$phosphopeptide.heavy.ratio)])
        df$SRM <- (df$value - median(x = df$value, na.rm = T))/IQR(x = df$value, na.rm = T)
        cap <- max(abs(df$SRM), na.rm = T)
        p3 = ggplot(df)
        p3 = p3 + geom_tile(aes(x=sample, y= paste0(sub, ":", rsd, "(SRM)"), fill= SRM, color = outlier_type, width=0.7, height=0.7), size=0.5)
        p3 = p3 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
        p3 = p3 + scale_fill_gradientn(name= "SRM phosphopeptide.heavy.ratio", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
        p3 = p3 + theme_bw() + theme_nogrid()
        p3 = p3 + theme(axis.text.x = element_blank())
        p3 = p3 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                        strip.background = element_blank(), strip.text.x = element_blank())
        p3 = p3 + theme(axis.ticks=element_blank(), 
                        legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
        p3 = p3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
        p3 = p3 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
        
        gb0 <- ggplot_build(p0)
        g0 <- ggplot_gtable(gb0)
        n0 <- length(gb0$layout$panel_ranges[[1]]$y.labels)
        
        gb1 <- ggplot_build(p1)
        gA <- ggplot_gtable(gb1)
        g <- gtable:::rbind_gtable(g0, gA, "last")
        n1 <- length(gb1$layout$panel_ranges[[1]]$y.labels)
        
        
        gb2 <- ggplot_build(p2)
        gB <- ggplot_gtable(gb2)
        g <- gtable:::rbind_gtable(g, gB, "last")
        n2 <- length(gb2$layout$panel_ranges[[1]]$y.labels)
        
        
        gb3 <- ggplot_build(p3)
        gC <- ggplot_gtable(gb3)
        g <- gtable:::rbind_gtable(g, gC, "last")
        n3 <- length(gb3$layout$panel_ranges[[1]]$y.labels)
        
        panels <- g$layout$t[grep("panel", g$layout$name)]
        g$layout$name[grep("panel", g$layout$name)]
        
        ## adjust height of each panel
        g$heights[panels] <- unit(x = c(n0, n0, n1, n1, n2, n2, n3, n3), units = "null")
        fn = paste0(makeOutDir(), cancer, "_", kinase, "_", sub, "_", rsd ,"_SRM.pdf")
        grid.newpage()
        pdf(fn, height = 1.5, width = 15)
        grid.draw(g)
        dev.off()
    }
    }
  }
}


# CO just SRM samples --------------------------------------------------------------------
colnames(trans_tab)[5] <- "peptide"
common_columns <- intersect(colnames(trans_tab), colnames(cis_tab))
tab <- rbind(trans_tab[, common_columns], cis_tab[, common_columns])
for (cancer in c("CO")) {
  pp_cancer <- pp_tab[, c("Peptide.Sequence", "Protein.Name", colnames(pp_tab)[grepl(pattern = cancer, x = colnames(pp_tab))])]
  tab_cancer <- tab[tab$cancer_type ==  cancer,]
  samples <- unique(tab_cancer[, c("Specimen_Label", "Participant_ID")])
  samples <- merge(samples, clinical[, c("Specimen.Label", "Participant.ID", "tumor_normal")], by.x = c("Specimen_Label"), by.y = c("Specimen.Label"), all.x = T)
  samples_dup <- unique(samples$Participant_ID[duplicated(samples$Participant_ID)])
  samples_nondup <- samples[!(samples %in% samples_dup),]
  
  ## input global data
  phog <- fread(input = paste(inputD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged.txt",sep=""), data.table = F)
  
  pho <- fread(input = paste(inputD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_replicate_averaged.txt",sep=""), data.table = F)
  pho <- cbind(formatPhosphosite(phosphosite_vector = pho$Phosphosite, pho$Gene), pho[,!(colnames(pho) %in% c("Gene", "Phosphosite"))])
  
  for (peptide in unique(pp_cancer$Peptide.Sequence)) {
    pp_cancer_pep <- pp_cancer[pp_cancer$Peptide.Sequence == toupper(peptide),]
    tab_cancer_pep <- tab_cancer[grepl(pattern = peptide, x = tab_cancer$peptide, ignore.case = T),]
    pp_cancer_pep.m <- melt(pp_cancer_pep, id.vars = c("Peptide.Sequence", "Protein.Name"))
    pp_cancer_pep.m$Participant_ID <- str_split_fixed(string = pp_cancer_pep.m$variable, pattern = "X", 2)[,2]
    pp_cancer_pep.m.nodup <- pp_cancer_pep.m[!(pp_cancer_pep.m$Participant_ID %in% samples_dup),]
    pp_cancer_pep.m.nodup <- merge(pp_cancer_pep.m.nodup, samples, by = c("Participant_ID"), all.x = T)
    
    sub <- unique(tab_cancer_pep$SUBSTRATE)
    rsd <- unique(tab_cancer_pep$SUB_MOD_RSD)
    pho_sub <- pho[pho$SUBSTRATE == sub & pho$SUB_MOD_RSD == rsd,]
    if (nrow(pho_sub) > 0 ){
      pho_sub.m <- melt(pho_sub, id.vars = c("SUBSTRATE", "transcript", "SUB_MOD_RSD"))
      colnames(pho_sub.m) <- c("substrate", "transcript", "SUB_MOD_RSD", "sample", "pho_sub")
      kinases <- unique(tab_cancer_pep$KINASE)
      for (kinase in kinases) {
        phog_kin <- phog[phog$Gene == kinase,]
        phog_kin.m <- melt(phog_kin, id.vars = c("Gene"))
        colnames(phog_kin.m) <- c("kinase", "sample", "phog_kin")
        df <- merge(pho_sub.m, phog_kin.m, all = T)
        df <- merge(df, pp_cancer_pep.m.nodup, by.x = c("sample"), by.y = c("Specimen_Label"), all = T)
        
        df$tumor_normal <- NULL
        df <- merge(df, clinical[, c("Specimen.Label", "tumor_normal")], by.x = c("sample"), by.y = c("Specimen.Label"), all.x = T)
        df <- df[!is.na(df$tumor_normal),]
        df$kinase <- kinase
        
        df <- merge(df, tab_cancer_pep, by.x = c("sample"), by.y = c("Specimen_Label"), all.x = T)
        df$value[df$value == "#N/A"] <- NA
        df$value <- as.numeric(as.vector(df$value))
        df$value[is.na(df$value) & !is.na(df$phosphopeptide.heavy.ratio)] <- as.numeric(df$phosphopeptide.heavy.ratio[is.na(df$value) & !is.na(df$phosphopeptide.heavy.ratio)])
        df$SRM <- (df$value - median(x = df$value, na.rm = T))/IQR(x = df$value, na.rm = T)
        
        df <- df[!is.na(df$SRM),]
        if (nrow(df) > 0 & length(unique(df$sample)) > 1) {
          df$sample <- as.vector(df$sample)
          df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$phog_kin)])
          df$tumor_normal <- factor(df$tumor_normal, levels = c("Tumor", "Normal"))
          
          p0 = ggplot(df)
          p0 = p0 + geom_tile(aes(x=sample, y= "Tumor (purple); Normal (green)", fill= tumor_normal, color = outlier_type, width=0.7, height=0.7), size=0.5)
          p0 = p0 + scale_fill_manual(values = c("Tumor" = "#6A3D9A", "Normal" = "#33A02C"))
          p0 = p0 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
          p0 = p0 + theme_bw() + theme_nogrid()
          p0 = p0 + theme(axis.text.x = element_blank())
          p0 = p0 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                          strip.background = element_blank())
          p0 = p0 + theme(axis.ticks=element_blank(), 
                          legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
          p0 = p0 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
          p0 = p0 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
          
          cap <- min(max(abs(df$phog_kin), na.rm = T),2)
          df$phog_kin_capped <- df$phog_kin
          df$phog_kin_capped[df$phog_kin > cap] <- cap
          df$phog_kin_capped[df$phog_kin < (-cap)] <- (-cap)
          p1 = ggplot(df)
          p1 = p1 + geom_tile(aes(x=sample, y= paste0(kinase, " phospho (global)"), fill= phog_kin_capped, color = outlier_type, width=0.7, height=0.7), size=0.5)
          p1 = p1 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
          p1 = p1 + scale_fill_gradientn(name= "global kinase phosphorylation abundance", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
          p1 = p1 + theme_bw() + theme_nogrid()
          p1 = p1 + theme(axis.text.x = element_blank())
          p1 = p1 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                          strip.background = element_blank())
          p1 = p1 + theme(axis.ticks=element_blank(), 
                          legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
          p1 = p1 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
          p1 = p1 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
          
          cap <- min(max(abs(df$pho_sub), na.rm = T),2)
          df$pho_sub_capped <- df$pho_sub
          df$pho_sub_capped[df$pho_sub > cap] <- cap
          df$pho_sub_capped[df$pho_sub < (-cap)] <- (-cap)
          p2 = ggplot(df)
          p2 = p2 + geom_tile(aes(x=sample, y= paste0(sub, ":", rsd, "(global)"), fill= pho_sub_capped, color = outlier_type, width=0.7, height=0.7), size=0.5)
          p2 = p2 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
          p2 = p2 + scale_fill_gradientn(name= "global substrate phosphosite abundance", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
          p2 = p2 + theme_bw() + theme_nogrid()
          p2 = p2 + theme(axis.text.x = element_blank())
          p2 = p2 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                          strip.background = element_blank(), strip.text.x = element_blank())
          p2 = p2 + theme(axis.ticks=element_blank(), 
                          legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
          p2 = p2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
          p2 = p2 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
          
          cap <- max(abs(df$SRM), na.rm = T)
          p3 = ggplot(df)
          p3 = p3 + geom_tile(aes(x=sample, y= paste0(sub, ":", rsd, "(SRM)"), fill= SRM, color = outlier_type, width=0.7, height=0.7), size=0.5)
          p3 = p3 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
          p3 = p3 + scale_fill_gradientn(name= "SRM phosphopeptide.heavy.ratio", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
          p3 = p3 + theme_bw() + theme_nogrid()
          p3 = p3 + theme(axis.text.x = element_blank())
          p3 = p3 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                          strip.background = element_blank(), strip.text.x = element_blank())
          p3 = p3 + theme(axis.ticks=element_blank(), 
                          legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
          p3 = p3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
          p3 = p3 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
          
          gb0 <- ggplot_build(p0)
          g0 <- ggplot_gtable(gb0)
          n0 <- length(gb0$layout$panel_ranges[[1]]$y.labels)
          
          gb1 <- ggplot_build(p1)
          gA <- ggplot_gtable(gb1)
          g <- gtable:::rbind_gtable(g0, gA, "last")
          n1 <- length(gb1$layout$panel_ranges[[1]]$y.labels)
          
          
          gb2 <- ggplot_build(p2)
          gB <- ggplot_gtable(gb2)
          g <- gtable:::rbind_gtable(g, gB, "last")
          n2 <- length(gb2$layout$panel_ranges[[1]]$y.labels)
          
          
          gb3 <- ggplot_build(p3)
          gC <- ggplot_gtable(gb3)
          g <- gtable:::rbind_gtable(g, gC, "last")
          n3 <- length(gb3$layout$panel_ranges[[1]]$y.labels)
          
          panels <- g$layout$t[grep("panel", g$layout$name)]
          g$layout$name[grep("panel", g$layout$name)]
          
          ## adjust height of each panel
          g$heights[panels] <- unit(x = c(n0, n0, n1, n1, n2, n2, n3, n3), units = "null")
          subdir <- paste0(makeOutDir(), cancer, "_justSRM_samples/")
          dir.create(subdir)
          fn = paste0(subdir, cancer, "_", kinase, "_", sub, "_", rsd ,"_justSRMsamples.pdf")
          grid.newpage()
          pdf(fn, height = 1.5, width = 6)
          grid.draw(g)
          dev.off()
        }
      }
    }
  }
}


# OV just SRM samples --------------------------------------------------------------------
colnames(trans_tab)[5] <- "peptide"
common_columns <- intersect(colnames(trans_tab), colnames(cis_tab))
tab <- rbind(trans_tab[, common_columns], cis_tab[, common_columns])
for (cancer in c("OV")) {
  pp_cancer <- pp_tab[, c("Peptide.Sequence", "Protein.Name", colnames(pp_tab)[grepl(pattern = cancer, x = colnames(pp_tab))])]
  tab_cancer <- tab[tab$cancer_type ==  cancer,]
  samples <- unique(tab_cancer[, c("Specimen_Label", "Participant_ID")])
  samples <- merge(samples, clinical[, c("Specimen.Label", "Participant.ID", "tumor_normal")], by.x = c("Specimen_Label"), by.y = c("Specimen.Label"), all.x = T)
  samples_dup <- unique(samples$Participant_ID[duplicated(samples$Participant_ID)])
  samples_nondup <- samples[!(samples %in% samples_dup),]
  
  ## input global data
  phog <- fread(input = paste(inputD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged.txt",sep=""), data.table = F)
  
  pho <- fread(input = paste(inputD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_replicate_averaged.txt",sep=""), data.table = F)
  pho <- cbind(formatPhosphosite(phosphosite_vector = pho$Phosphosite, pho$Gene), pho[,!(colnames(pho) %in% c("Gene", "Phosphosite"))])
  
  for (peptide in unique(pp_cancer$Peptide.Sequence)) {
    pp_cancer_pep <- pp_cancer[pp_cancer$Peptide.Sequence == toupper(peptide),]
    tab_cancer_pep <- tab_cancer[grepl(pattern = peptide, x = tab_cancer$peptide, ignore.case = T),]
    pp_cancer_pep.m <- melt(pp_cancer_pep, id.vars = c("Peptide.Sequence", "Protein.Name"))
    pp_cancer_pep.m$Participant_ID <- str_split_fixed(string = pp_cancer_pep.m$variable, pattern = "X", 2)[,2]
    pp_cancer_pep.m.nodup <- pp_cancer_pep.m[!(pp_cancer_pep.m$Participant_ID %in% samples_dup),]
    pp_cancer_pep.m.nodup <- merge(pp_cancer_pep.m.nodup, samples, by = c("Participant_ID"), all.x = T)
    
    sub <- unique(tab_cancer_pep$SUBSTRATE)
    rsd <- unique(tab_cancer_pep$SUB_MOD_RSD)
    pho_sub <- pho[pho$SUBSTRATE == sub & pho$SUB_MOD_RSD == rsd,]
    if (nrow(pho_sub) > 0 ){
      pho_sub.m <- melt(pho_sub, id.vars = c("SUBSTRATE", "transcript", "SUB_MOD_RSD"))
      colnames(pho_sub.m) <- c("substrate", "transcript", "SUB_MOD_RSD", "sample", "pho_sub")
      kinases <- unique(tab_cancer_pep$KINASE)
      for (kinase in kinases) {
        phog_kin <- phog[phog$Gene == kinase,]
        phog_kin.m <- melt(phog_kin, id.vars = c("Gene"))
        colnames(phog_kin.m) <- c("kinase", "sample", "phog_kin")
        df <- merge(pho_sub.m, phog_kin.m, all = T)
        df <- merge(df, pp_cancer_pep.m.nodup, by.x = c("sample"), by.y = c("Specimen_Label"), all = T)
        
        df$tumor_normal <- NULL
        df <- merge(df, clinical[, c("Specimen.Label", "tumor_normal")], by.x = c("sample"), by.y = c("Specimen.Label"), all.x = T)
        df <- df[!is.na(df$tumor_normal),]
        df$kinase <- kinase
        
        df <- merge(df, tab_cancer_pep, by.x = c("sample"), by.y = c("Specimen_Label"), all.x = T)
        df$value[df$value == "#N/A"] <- NA
        df$value <- as.numeric(as.vector(df$value))
        df$value[is.na(df$value) & !is.na(df$phosphopeptide.heavy.ratio)] <- as.numeric(df$phosphopeptide.heavy.ratio[is.na(df$value) & !is.na(df$phosphopeptide.heavy.ratio)])
        df$SRM <- (df$value - median(x = df$value, na.rm = T))/IQR(x = df$value, na.rm = T)
        
        df <- df[!is.na(df$SRM),]
        if (nrow(df) > 0 & length(unique(df$sample)) > 1) {
          df$sample <- as.vector(df$sample)
          df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$phog_kin)])
          df$tumor_normal <- factor(df$tumor_normal, levels = c("Tumor", "Normal"))
          
          p0 = ggplot(df)
          p0 = p0 + geom_tile(aes(x=sample, y= "Tumor (purple); Normal (green)", fill= tumor_normal, color = outlier_type, width=0.7, height=0.7), size=0.5)
          p0 = p0 + scale_fill_manual(values = c("Tumor" = "#6A3D9A", "Normal" = "#33A02C"))
          p0 = p0 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
          p0 = p0 + theme_bw() + theme_nogrid()
          p0 = p0 + theme(axis.text.x = element_blank())
          p0 = p0 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                          strip.background = element_blank())
          p0 = p0 + theme(axis.ticks=element_blank(), 
                          legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
          p0 = p0 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
          p0 = p0 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
          
          cap <- min(max(abs(df$phog_kin), na.rm = T),2)
          df$phog_kin_capped <- df$phog_kin
          df$phog_kin_capped[df$phog_kin > cap] <- cap
          df$phog_kin_capped[df$phog_kin < (-cap)] <- (-cap)
          p1 = ggplot(df)
          p1 = p1 + geom_tile(aes(x=sample, y= paste0(kinase, " phospho (global)"), fill= phog_kin_capped, color = outlier_type, width=0.7, height=0.7), size=0.5)
          p1 = p1 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
          p1 = p1 + scale_fill_gradientn(name= "global kinase phosphorylation abundance", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
          p1 = p1 + theme_bw() + theme_nogrid()
          p1 = p1 + theme(axis.text.x = element_blank())
          p1 = p1 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                          strip.background = element_blank())
          p1 = p1 + theme(axis.ticks=element_blank(), 
                          legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
          p1 = p1 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
          p1 = p1 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
          
          cap <- min(max(abs(df$pho_sub), na.rm = T),2)
          df$pho_sub_capped <- df$pho_sub
          df$pho_sub_capped[df$pho_sub > cap] <- cap
          df$pho_sub_capped[df$pho_sub < (-cap)] <- (-cap)
          p2 = ggplot(df)
          p2 = p2 + geom_tile(aes(x=sample, y= paste0(sub, ":", rsd, "(global)"), fill= pho_sub_capped, color = outlier_type, width=0.7, height=0.7), size=0.5)
          p2 = p2 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
          p2 = p2 + scale_fill_gradientn(name= "global substrate phosphosite abundance", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
          p2 = p2 + theme_bw() + theme_nogrid()
          p2 = p2 + theme(axis.text.x = element_blank())
          p2 = p2 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                          strip.background = element_blank(), strip.text.x = element_blank())
          p2 = p2 + theme(axis.ticks=element_blank(), 
                          legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
          p2 = p2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
          p2 = p2 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
          
          cap <- max(abs(df$SRM), na.rm = T)
          p3 = ggplot(df)
          p3 = p3 + geom_tile(aes(x=sample, y= paste0(sub, ":", rsd, "(SRM)"), fill= SRM, color = outlier_type, width=0.7, height=0.7), size=0.5)
          p3 = p3 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
          p3 = p3 + scale_fill_gradientn(name= "SRM phosphopeptide.heavy.ratio", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
          p3 = p3 + theme_bw() + theme_nogrid()
          p3 = p3 + theme(axis.text.x = element_blank())
          p3 = p3 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                          strip.background = element_blank(), strip.text.x = element_blank())
          p3 = p3 + theme(axis.ticks=element_blank(), 
                          legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
          p3 = p3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
          p3 = p3 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
          
          gb0 <- ggplot_build(p0)
          g0 <- ggplot_gtable(gb0)
          n0 <- length(gb0$layout$panel_ranges[[1]]$y.labels)
          
          gb1 <- ggplot_build(p1)
          gA <- ggplot_gtable(gb1)
          g <- gtable:::rbind_gtable(g0, gA, "last")
          n1 <- length(gb1$layout$panel_ranges[[1]]$y.labels)
          
          
          gb2 <- ggplot_build(p2)
          gB <- ggplot_gtable(gb2)
          g <- gtable:::rbind_gtable(g, gB, "last")
          n2 <- length(gb2$layout$panel_ranges[[1]]$y.labels)
          
          
          gb3 <- ggplot_build(p3)
          gC <- ggplot_gtable(gb3)
          g <- gtable:::rbind_gtable(g, gC, "last")
          n3 <- length(gb3$layout$panel_ranges[[1]]$y.labels)
          
          panels <- g$layout$t[grep("panel", g$layout$name)]
          g$layout$name[grep("panel", g$layout$name)]
          
          ## adjust height of each panel
          g$heights[panels] <- unit(x = c(n0, n0, n1, n1, n2, n2, n3, n3), units = "null")
          subdir <- paste0(makeOutDir(), cancer, "_justSRM_samples/")
          dir.create(subdir)
          fn = paste0(subdir, cancer, "_", kinase, "_", sub, "_", rsd ,"_justSRMsamples.pdf")
          grid.newpage()
          pdf(fn, height = 1.5, width = 6)
          grid.draw(g)
          dev.off()
        }
      }
    }
  }
}

# BRCA just SRM samples --------------------------------------------------------------------
colnames(trans_tab)[5] <- "peptide"
common_columns <- intersect(colnames(trans_tab), colnames(cis_tab))
tab <- rbind(trans_tab[, common_columns], cis_tab[, common_columns])
for (cancer in c("BRCA")) {
  pp_cancer <- pp_tab[, c("Peptide.Sequence", "Protein.Name", colnames(pp_tab)[grepl(pattern = "BR", x = colnames(pp_tab))])]
  tab_cancer <- tab[tab$cancer_type ==  cancer,]
  samples <- unique(tab_cancer[, c("Specimen_Label", "Participant_ID")])
  samples <- merge(samples, clinical[, c("Specimen.Label", "Participant.ID", "tumor_normal")], by.x = c("Specimen_Label"), by.y = c("Specimen.Label"), all.x = T)
  samples_dup <- unique(samples$Participant_ID[duplicated(samples$Participant_ID)])
  samples_nondup <- samples[!(samples %in% samples_dup),]
  
  ## input global data
  phog <- fread(input = paste(inputD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged.txt",sep=""), data.table = F)
  
  pho <- fread(input = paste(inputD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_replicate_averaged.txt",sep=""), data.table = F)
  pho <- cbind(formatPhosphosite(phosphosite_vector = pho$Phosphosite, pho$Gene), pho[,!(colnames(pho) %in% c("Gene", "Phosphosite"))])
  
  for (peptide in unique(pp_cancer$Peptide.Sequence)) {
    pp_cancer_pep <- pp_cancer[pp_cancer$Peptide.Sequence == toupper(peptide),]
    tab_cancer_pep <- tab_cancer[grepl(pattern = peptide, x = tab_cancer$peptide, ignore.case = T),]
    pp_cancer_pep.m <- melt(pp_cancer_pep, id.vars = c("Peptide.Sequence", "Protein.Name"))
    pp_cancer_pep.m$Participant_ID <- str_split_fixed(string = pp_cancer_pep.m$variable, pattern = "X", 2)[,2]
    pp_cancer_pep.m.nodup <- pp_cancer_pep.m[!(pp_cancer_pep.m$Participant_ID %in% samples_dup),]
    pp_cancer_pep.m.nodup <- merge(pp_cancer_pep.m.nodup, samples, by = c("Participant_ID"), all.x = T)
    
    sub <- unique(tab_cancer_pep$SUBSTRATE)
    rsd <- unique(tab_cancer_pep$SUB_MOD_RSD)
    pho_sub <- pho[pho$SUBSTRATE == sub & pho$SUB_MOD_RSD == rsd,]
    if (nrow(pho_sub) > 0 ){
      pho_sub.m <- melt(pho_sub, id.vars = c("SUBSTRATE", "transcript", "SUB_MOD_RSD"))
      colnames(pho_sub.m) <- c("substrate", "transcript", "SUB_MOD_RSD", "sample", "pho_sub")
      kinases <- unique(tab_cancer_pep$KINASE)
      for (kinase in kinases) {
        phog_kin <- phog[phog$Gene == kinase,]
        phog_kin.m <- melt(phog_kin, id.vars = c("Gene"))
        colnames(phog_kin.m) <- c("kinase", "sample", "phog_kin")
        df <- merge(pho_sub.m, phog_kin.m, all = T)
        df <- merge(df, pp_cancer_pep.m.nodup, by.x = c("sample"), by.y = c("Specimen_Label"), all = T)
        
        df$tumor_normal <- NULL
        df <- merge(df, clinical[, c("Specimen.Label", "tumor_normal")], by.x = c("sample"), by.y = c("Specimen.Label"), all.x = T)
        df <- df[!is.na(df$tumor_normal),]
        df$kinase <- kinase
        
        df <- merge(df, tab_cancer_pep, by.x = c("sample"), by.y = c("Specimen_Label"), all.x = T)
        df$value[df$value == "#N/A"] <- NA
        df$value <- as.numeric(as.vector(df$value))
        df$value[is.na(df$value) & !is.na(df$phosphopeptide.heavy.ratio)] <- as.numeric(df$phosphopeptide.heavy.ratio[is.na(df$value) & !is.na(df$phosphopeptide.heavy.ratio)])
        df$SRM <- (df$value - median(x = df$value, na.rm = T))/IQR(x = df$value, na.rm = T)
        
        df <- df[!is.na(df$SRM),]
        if (nrow(df) > 0 & length(unique(df$sample)) > 1) {
          df$sample <- as.vector(df$sample)
          df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$phog_kin)])
          df$tumor_normal <- factor(df$tumor_normal, levels = c("Tumor", "Normal"))
          
          p0 = ggplot(df)
          p0 = p0 + geom_tile(aes(x=sample, y= "Tumor (purple); Normal (green)", fill= tumor_normal, color = outlier_type, width=0.7, height=0.7), size=0.5)
          p0 = p0 + scale_fill_manual(values = c("Tumor" = "#6A3D9A", "Normal" = "#33A02C"))
          p0 = p0 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
          p0 = p0 + theme_bw() + theme_nogrid()
          p0 = p0 + theme(axis.text.x = element_blank())
          p0 = p0 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                          strip.background = element_blank())
          p0 = p0 + theme(axis.ticks=element_blank(), 
                          legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
          p0 = p0 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
          p0 = p0 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
          
          cap <- min(max(abs(df$phog_kin), na.rm = T),2)
          df$phog_kin_capped <- df$phog_kin
          df$phog_kin_capped[df$phog_kin > cap] <- cap
          df$phog_kin_capped[df$phog_kin < (-cap)] <- (-cap)
          p1 = ggplot(df)
          p1 = p1 + geom_tile(aes(x=sample, y= paste0(kinase, " phospho (global)"), fill= phog_kin_capped, color = outlier_type, width=0.7, height=0.7), size=0.5)
          p1 = p1 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
          p1 = p1 + scale_fill_gradientn(name= "global kinase phosphorylation abundance", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
          p1 = p1 + theme_bw() + theme_nogrid()
          p1 = p1 + theme(axis.text.x = element_blank())
          p1 = p1 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                          strip.background = element_blank())
          p1 = p1 + theme(axis.ticks=element_blank(), 
                          legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
          p1 = p1 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
          p1 = p1 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
          
          cap <- min(max(abs(df$pho_sub), na.rm = T),2)
          df$pho_sub_capped <- df$pho_sub
          df$pho_sub_capped[df$pho_sub > cap] <- cap
          df$pho_sub_capped[df$pho_sub < (-cap)] <- (-cap)
          p2 = ggplot(df)
          p2 = p2 + geom_tile(aes(x=sample, y= paste0(sub, ":", rsd, "(global)"), fill= pho_sub_capped, color = outlier_type, width=0.7, height=0.7), size=0.5)
          p2 = p2 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
          p2 = p2 + scale_fill_gradientn(name= "global substrate phosphosite abundance", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
          p2 = p2 + theme_bw() + theme_nogrid()
          p2 = p2 + theme(axis.text.x = element_blank())
          p2 = p2 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                          strip.background = element_blank(), strip.text.x = element_blank())
          p2 = p2 + theme(axis.ticks=element_blank(), 
                          legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
          p2 = p2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
          p2 = p2 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
          
          cap <- max(abs(df$SRM), na.rm = T)
          p3 = ggplot(df)
          p3 = p3 + geom_tile(aes(x=sample, y= paste0(sub, ":", rsd, "(SRM)"), fill= SRM, color = outlier_type, width=0.7, height=0.7), size=0.5)
          p3 = p3 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
          p3 = p3 + scale_fill_gradientn(name= "SRM phosphopeptide.heavy.ratio", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
          p3 = p3 + theme_bw() + theme_nogrid()
          p3 = p3 + theme(axis.text.x = element_blank())
          p3 = p3 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7),
                          strip.background = element_blank(), strip.text.x = element_blank())
          p3 = p3 + theme(axis.ticks=element_blank(), 
                          legend.text = element_text(size = 5), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
          p3 = p3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
          p3 = p3 +  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1)))
          
          gb0 <- ggplot_build(p0)
          g0 <- ggplot_gtable(gb0)
          n0 <- length(gb0$layout$panel_ranges[[1]]$y.labels)
          
          gb1 <- ggplot_build(p1)
          gA <- ggplot_gtable(gb1)
          g <- gtable:::rbind_gtable(g0, gA, "last")
          n1 <- length(gb1$layout$panel_ranges[[1]]$y.labels)
          
          
          gb2 <- ggplot_build(p2)
          gB <- ggplot_gtable(gb2)
          g <- gtable:::rbind_gtable(g, gB, "last")
          n2 <- length(gb2$layout$panel_ranges[[1]]$y.labels)
          
          
          gb3 <- ggplot_build(p3)
          gC <- ggplot_gtable(gb3)
          g <- gtable:::rbind_gtable(g, gC, "last")
          n3 <- length(gb3$layout$panel_ranges[[1]]$y.labels)
          
          panels <- g$layout$t[grep("panel", g$layout$name)]
          g$layout$name[grep("panel", g$layout$name)]
          
          ## adjust height of each panel
          g$heights[panels] <- unit(x = c(n0, n0, n1, n1, n2, n2, n3, n3), units = "null")
          subdir <- paste0(makeOutDir(), cancer, "_justSRM_samples/")
          dir.create(subdir)
          fn = paste0(subdir, cancer, "_", kinase, "_", sub, "_", rsd ,"_justSRMsamples.pdf")
          grid.newpage()
          pdf(fn, height = 1.5, width = 6)
          grid.draw(g)
          dev.off()
        }
      }
    }
  }
}
