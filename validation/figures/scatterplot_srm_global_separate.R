# Yige Wu @ WashU 2018 Apr
## make heatmap showing the ks outlier samples really agree regarding the order of abundance


# source ------------------------------------------------------------------
library(readxl)
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')

library(GGally)

# inputs ------------------------------------------------------------------
## input the date of receiving the data
email_date <- "2018-04-20"

## input the sample ID mapping info
clinical <- loadSampMap()

## input the SRM data
# cis_tab <- read_excel("~/Box Sync/cptac2p/tao/2018-03-29/Druggable-kinase-outlier-CPTAC2-SRM_03292018_omnipath_result_commented.xlsx",
#                       sheet = "cis")
# trans_tab <- read_excel("~/Box Sync/cptac2p/tao/2018-03-29/Druggable-kinase-outlier-CPTAC2-SRM_03292018_omnipath_result_commented.xlsx",
#                       sheet = "trans")
# pp_tab <- read_excel(path = "~/Box Sync/cptac2p/tao/2018-03-29/Phosphopeptide_Crosstab.xlsx")

cis_tab <- read_excel("~/Box Sync/cptac2p/tao/2018-04-20/Druggable-kinase-outlier-CPTAC2prospective-SRM_04192018.xlsx",
                      sheet = "cis")
trans_tab <- read_excel("~/Box Sync/cptac2p/tao/2018-04-20/Druggable-kinase-outlier-CPTAC2prospective-SRM_04192018.xlsx",
                        sheet = "trans")
pp_tab <- read_excel(path = "~/Box Sync/cptac2p/tao/2018-04-20/Phosphopeptide_04192018_Crosstab.xlsx")
up_tab <- read_excel(path = "~/Box Sync/cptac2p/tao/2018-04-20/Unmodified peptide_04192018_Crosstab.xlsx")

cis_tab <- data.frame(cis_tab)
trans_tab <- data.frame(trans_tab)
pp_tab <- data.frame(pp_tab)
up_tab <- data.frame(up_tab)

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

# trans --------------------------------------------------------------------
colnames(trans_tab)[5] <- "peptide"
common_columns <- intersect(colnames(trans_tab), colnames(cis_tab))
tab <- rbind(trans_tab[, common_columns], cis_tab[, common_columns])
resultDnow <- makeOutDir(resultD = resultD)
subdir <- paste0(resultDnow, "trans", "/")
dir.create(subdir)
for (cancer in cancers_sort) {
# for (cancer in c("BRCA")) {
  pp_cancer <- pp_tab[, c("Peptide.Sequence", "Protein.Name", colnames(pp_tab)[grepl(pattern = substr(cancer, start = 1, stop = 2), x = colnames(pp_tab))])]
  tab_cancer <- tab[(as.vector(tab$KINASE) != as.vector(tab$SUBSTRATE)),]
  
  ## input global data
  phog <- loadPhosphoproteinNormalizedTumor(cancer = cancer)
  pho <- loadPhosphositeNormalizedTumor(cancer = cancer)
  pho <- cbind(formatPhosphosite(phosphosite_vector = pho$Phosphosite, pho$Gene), pho[,!(colnames(pho) %in% c("Gene", "Phosphosite"))])
  
  for (peptide in unique(pp_cancer$Peptide.Sequence)) {
    pp_cancer_pep <- pp_cancer[pp_cancer$Peptide.Sequence == toupper(peptide),]
    tab_cancer_pep <- tab_cancer[grepl(pattern = peptide, x = tab_cancer$peptide, ignore.case = T),]
    tab_pep <- tab[grepl(pattern = peptide, x = tab$peptide, ignore.case = T),]
    
    pp_cancer_pep.m <- melt(pp_cancer_pep, id.vars = c("Peptide.Sequence", "Protein.Name"))
    pp_cancer_pep.m$Participant_ID <- str_split_fixed(string = pp_cancer_pep.m$variable, pattern = "X", 2)[,2]

    sub <- unique(tab_cancer_pep$SUBSTRATE)
    rsd <- unique(tab_cancer_pep$SUB_MOD_RSD)
    pho_sub <- pho[pho$SUBSTRATE == sub & pho$SUB_MOD_RSD == rsd,]
    if (nrow(pho_sub) > 0 ){
      pho_sub.m <- melt(pho_sub, id.vars = c("SUBSTRATE", "transcript", "SUB_MOD_RSD"))
      colnames(pho_sub.m) <- c("substrate", "transcript", "SUB_MOD_RSD", "sample", "pho_sub")
      kinases <- unique(tab_pep$KINASE)
      kinases <- kinases[kinases!=sub]
      for (kinase in kinases) {
        subdir1 <- paste0(subdir, kinase, "_", sub, "_", rsd, "/")
        dir.create(subdir1)
        
        subdir2 <- paste0(subdir1, "phospho-peptide_",  peptide, "/")
        dir.create(subdir2)
        
        phog_kin <- phog[phog$Gene == kinase,]
        phog_kin.m <- melt(phog_kin, id.vars = c("Gene"))
        colnames(phog_kin.m) <- c("kinase", "sample", "phog_kin")
        phog_kin.m <- merge(phog_kin.m, clinical[, c("Specimen.Label", "Participant.ID", "tumor_normal")], by.x = c("sample"), by.y = c("Specimen.Label"), all.x = T)
        
        df1 <- merge(pho_sub.m[,c("sample", "pho_sub")], phog_kin.m[,c("sample", "phog_kin")], all = T)
        df1$technology <- "global"
        df1$Participant_ID <- sampID2partID(sampleID_vector = as.vector(df1$sample), sample_map = clinical)
        
        df2 <- merge(pp_cancer_pep.m[,c("Participant_ID", "value")], phog_kin.m[phog_kin.m$tumor_normal == "Tumor",c("Participant.ID", "phog_kin")], by.y = c("Participant.ID"), by.x = c("Participant_ID"), all = T)
        df2$value[df2$value == "#N/A"] <- NA
        df2$value <- as.numeric(as.vector(df2$value))
        df2$value[is.na(df2$value) & !is.na(df2$phosphopeptide.heavy.ratio)] <- as.numeric(df2$phosphopeptide.heavy.ratio[is.na(df2$value) & !is.na(df2$phosphopeptide.heavy.ratio)])
        colnames(df2) <- c("Participant_ID", "pho_sub", "phog_kin")
        df2$technology <- "SRM"
        df2$tumor_normal <- "Tumor"
        col_tmp <- intersect(colnames(df1), colnames(df2))
        df <- rbind(df1[,col_tmp], df2[,col_tmp])
        df <- df[!is.na(remove_outliers(df$phog_kin, out_thres = 1.5)),]

        ## reshape the data.frame
        df_cor <- df[df$technology == "global", c("phog_kin", "pho_sub")]
        colnames(df_cor)[ncol(df_cor)] <- "pho_sub.global"
        df_cor <- merge(df_cor, df[df$technology == "SRM", c("phog_kin", "pho_sub")], all.x = T)
        colnames(df_cor)[ncol(df_cor)] <- "pho_sub.srm"
        
        p = ggplot(df)
        p = p + geom_point(aes(x=phog_kin, y=pho_sub), alpha=0.9, stroke = 0)
        p <- p + geom_smooth(formula = y ~ x, mapping = aes(x = phog_kin, y = pho_sub), method = lm, alpha = 0.4)
        p <- p + facet_grid(technology~., scales = "free_y", space = "fixed")
        p = p + labs(x = paste0("global phosphorylaion abundance for ", kinase),
                     y=paste0("SRM/global phosphorylation abundance for ", sub, " ", rsd))
        p = p + theme(axis.title = element_text(size=10),
                      axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5),
                      axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
        p = p + theme(title = element_text(size = 8))
        p <- p + theme_bw()
        p
        fn = paste0(subdir2, cancer, "_", kinase, "_", sub, "_", rsd ,"_global~SRM.pdf")
        ggsave(file=fn, height=6, width=5, useDingbats=FALSE)
        
      }
    }
  }
}


# cis --------------------------------------------------------------------
colnames(trans_tab)[5] <- "peptide"
common_columns <- intersect(colnames(trans_tab), colnames(cis_tab))
tab <- rbind(trans_tab[, common_columns], cis_tab[, common_columns])
resultDnow <- makeOutDir(resultD = resultD)
subdir <- paste0(resultDnow, "cis", "/")
dir.create(subdir)
for (cancer in c("CO", "OV", "BRCA")) {
# for (cancer in c("OV")) {

  ## input cross tab for SRM values
  pp_cancer <- pp_tab[, c("Peptide.Sequence", "Protein.Name", colnames(pp_tab)[grepl(pattern = substr(cancer, start = 1, stop = 2), x = colnames(pp_tab))])]
  up_cancer <- up_tab[, c("Peptide.Sequence", "Protein.Name", colnames(pp_tab)[grepl(pattern = substr(cancer, start = 1, stop = 2), x = colnames(up_tab))])]
  
  tab_cancer <- tab[tab$cancer_type ==  cancer,]
  
  ## input global data
  phog <- loadPhosphoproteinNormalizedTumor(cancer = cancer)
  pho <- loadPhosphositeNormalizedTumor(cancer = cancer)
  pho <- cbind(formatPhosphosite(phosphosite_vector = pho$Phosphosite, pho$Gene), pho[,!(colnames(pho) %in% c("Gene", "Phosphosite"))])
  pro <- loadProteinNormalizedTumor(cancer = cancer)
  
  for (peptide in unique(pp_cancer$Peptide.Sequence)[4]) {
    pp_cancer_pep <- pp_cancer[pp_cancer$Peptide.Sequence == toupper(peptide),]
    pp_cancer_pep.m <- melt(pp_cancer_pep, id.vars = c("Peptide.Sequence", "Protein.Name"))
    pp_cancer_pep.m$Participant_ID <- str_split_fixed(string = pp_cancer_pep.m$variable, pattern = "X", 2)[,2]

    tab_cancer_pep <- tab_cancer[grepl(pattern = peptide, x = tab_cancer$peptide, ignore.case = T),]
    tab_pep <- tab[grepl(pattern = peptide, x = tab$peptide, ignore.case = T),]
    
    sub <- unique(tab_pep$SUBSTRATE)
    rsd <- unique(tab_pep$SUB_MOD_RSD)
    pho_sub <- pho[pho$SUBSTRATE == sub & pho$SUB_MOD_RSD == rsd,]
    
    if (nrow(pho_sub) > 0 ){
      pho_sub.m <- melt(pho_sub, id.vars = c("SUBSTRATE", "transcript", "SUB_MOD_RSD"))
      colnames(pho_sub.m) <- c("substrate", "transcript", "SUB_MOD_RSD", "sample", "pho_sub")
      if (sub %in% unique(tab_pep$KINASE)) {
        kinase <- sub
        subdir1 <- paste0(subdir, kinase, "_", sub, "_", rsd, "/")
        dir.create(subdir1)
        
        subdir2 <- paste0(subdir1, "phospho-peptide_",  peptide, "/")
        dir.create(subdir2)
        
        pro_kin <- pro[pro$Gene == kinase,]
        pro_kin.m <- melt(pro_kin, id.vars = c("Gene"))
        colnames(pro_kin.m) <- c("kinase", "sample", "pro_kin")
        pro_kin.m <- merge(pro_kin.m, clinical[, c("Specimen.Label", "Participant.ID", "tumor_normal")], by.x = c("sample"), by.y = c("Specimen.Label"), all.x = T)
        
        df1 <- merge(pho_sub.m[,c("sample", "pho_sub")], pro_kin.m[,c("sample", "pro_kin")], all = T)
        colnames(df1) <- c("sample", "pho_sub.global", "pro_kin.global")
        df1$Participant_ID <- sampID2partID(sampleID_vector = as.vector(df1$sample), sample_map = clinical)
        
        df2 <- merge(df1, pp_cancer_pep.m[,c("Participant_ID", "value")], by = c("Participant_ID"), all = T)
        df2$value[df2$value == "#N/A"] <- NA
        df2$value <- as.numeric(as.vector(df2$value))
        colnames(df2)[colnames(df2) == "value"] <- c("pho_sub.srm")
        
        ## try different unmodified peptide for protein level if possible
        for (upeptide in up_cancer$Peptide.Sequence[up_cancer$Protein.Name == pp_cancer_pep$Protein.Name]) {
          subdir3 <- paste0(subdir2, "unmodified-peptide_", upeptide, "/")
          dir.create(subdir3)
          
          up_cancer_pep <- up_cancer[up_cancer$Peptide.Sequence == upeptide,]
          up_cancer_pep.m <- melt(up_cancer_pep, id.vars = c("Peptide.Sequence", "Protein.Name"))
          up_cancer_pep.m$Participant_ID <- str_split_fixed(string = up_cancer_pep.m$variable, pattern = "X", 2)[,2]

          df3 <- merge(df2, up_cancer_pep.m[,c("Participant_ID", "value")], by = c("Participant_ID"), all = T)
          df3$value[df3$value == "#N/A"] <- NA
          df3$value <- as.numeric(as.vector(df3$value))
          colnames(df3)[colnames(df3) == "value"] <- c("pro_kin.srm")
          df3$srm_samples <- !is.na(df3$pho_sub.srm)

          
          # plot global correlation -------------------------------------------------
          df <- df3
          df$y <- remove_outliers(x = df$pho_sub.global, out_thres = 3, na.rm = T)
          df$y[!is.na(df$pho_sub.srm)] <- df$pho_sub.global[!is.na(df$pho_sub.srm)]
          df$x <- remove_outliers(x = df$pro_kin.global, out_thres = 3, na.rm = T)
          df$x[!is.na(df$pro_kin.srm)] <- df$pro_kin.global[!is.na(df$pro_kin.srm)]
          
          fn = paste0(subdir3, cancer, "_", kinase, "_", upeptide, "_", sub, "_", rsd, "_", peptide ,"_global.pdf")
          p <- ggplot()
          p <- p + geom_point(data = df, mapping = aes(x = x, y = y, color = srm_samples), alpha = 0.8)
          p <- p + scale_color_manual(values = c("TRUE" = set1[1], "FALSE" = "black"))
          p <- p + xlab(paste0(kinase, " protein abundance(global)")) + ylab(paste0(sub, " ", rsd, " abundance (global)"))
          p <- p + theme_bw()
          pdf(fn,  height=3, width=4.5, useDingbats = F)
          print(p)
          dev.off()
          
          df <- df3
          fn = paste0(subdir3, cancer, "_", kinase, "_", upeptide, "_", sub, "_", rsd, "_", peptide ,"_srm.pdf")
          p <- ggplot()
          p <- p + geom_point(data = df, mapping = aes(x = pro_kin.srm, y = pho_sub.srm, color = srm_samples), alpha = 0.8)
          p <- p + scale_color_manual(values = c("TRUE" = set1[1], "FALSE" = "black"))
          p <- p + xlab(paste0(kinase, " protein abundance(SRM)")) + ylab(paste0(sub, " ", rsd, " abundance (SRM)"))
          p <- p + theme_bw()
          pdf(fn,  height=3, width=4.5, useDingbats = F)
          print(p)
          dev.off()
          
        }
      }
    }
  }
}



