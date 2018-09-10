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
# trans_tab <- read_excel("~/Box Sync/cptac2p/tao/2018-03-29/Druggable-kinase-outlier-CPTAC2-SRM_03292018_omnipath_result_commented.xlsx",
#                       sheet = "trans")
# pp_tab <- read_excel(path = "~/Box Sync/cptac2p/tao/2018-03-29/Phosphopeptide_Crosstab.xlsx")

cis_tab <- read_excel("~/Box Sync/cptac2p/tao/2018-04-20/Druggable-kinase-outlier-CPTAC2prospective-SRM_04192018.xlsx",
                      sheet = "cis")
cis_tab <- data.frame(cis_tab)

trans_tab <- read_excel("~/Box Sync/cptac2p/tao/2018-04-20/Druggable-kinase-outlier-CPTAC2prospective-SRM_04192018.xlsx",
                        sheet = "trans")
trans_tab <- data.frame(trans_tab)

pp_tab <- read_excel(path = "~/Box Sync/cptac2p/tao/2018-04-20/Phosphopeptide_04192018_Crosstab.xlsx")
pp_tab <- data.frame(pp_tab)

up_tab <- read_excel(path = "~/Box Sync/cptac2p/tao/2018-04-20/Unmodified peptide_04192018_Crosstab.xlsx")
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
#for (cancer in c("CO", "OV", "BRCA")) {
for (cancer in c("BRCA")) {
    
  pp_cancer <- pp_tab[, c("Peptide.Sequence", "Protein.Name", colnames(pp_tab)[grepl(pattern = substr(cancer, start = 1, stop = 2), x = colnames(pp_tab))])]
  tab_cancer <- tab[tab$cancer_type ==  cancer,]
  
  ## input global data
  phog <- fread(input = paste(inputD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged.txt",sep=""), data.table = F)
  
  pho <- fread(input = paste(inputD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_replicate_averaged.txt",sep=""), data.table = F)
  pho <- cbind(formatPhosphosite(phosphosite_vector = pho$Phosphosite, pho$Gene), pho[,!(colnames(pho) %in% c("Gene", "Phosphosite"))])
  
  resultDnow <- makeOutDir()
  subdir <- paste0(resultDnow, cancer, "/")
  dir.create(subdir)
  for (peptide in unique(pp_cancer$Peptide.Sequence)) {
    pp_cancer_pep <- pp_cancer[pp_cancer$Peptide.Sequence == toupper(peptide),]
    tab_cancer_pep <- tab_cancer[grepl(pattern = peptide, x = tab_cancer$peptide, ignore.case = T),]
    pp_cancer_pep.m <- melt(pp_cancer_pep, id.vars = c("Peptide.Sequence", "Protein.Name"))
    pp_cancer_pep.m$Participant_ID <- str_split_fixed(string = pp_cancer_pep.m$variable, pattern = "X", 2)[,2]

    sub <- unique(tab_cancer_pep$SUBSTRATE)
    rsd <- unique(tab_cancer_pep$SUB_MOD_RSD)
    pho_sub <- pho[pho$SUBSTRATE == sub & pho$SUB_MOD_RSD == rsd,]
    if (nrow(pho_sub) > 0 ){
      pho_sub.m <- melt(pho_sub, id.vars = c("SUBSTRATE", "transcript", "SUB_MOD_RSD"))
      colnames(pho_sub.m) <- c("substrate", "transcript", "SUB_MOD_RSD", "sample", "pho_sub")
      kinases <- unique(tab_cancer_pep$KINASE)
      kinases <- kinases[kinases!=sub]
      for (kinase in kinases) {
        phog_kin <- phog[phog$Gene == kinase,]
        phog_kin.m <- melt(phog_kin, id.vars = c("Gene"))
        colnames(phog_kin.m) <- c("kinase", "sample", "phog_kin")
        phog_kin.m <- merge(phog_kin.m, clinical[, c("Specimen.Label", "Participant.ID", "tumor_normal")], by.x = c("sample"), by.y = c("Specimen.Label"), all.x = T)
        df1 <- merge(pho_sub.m[,c("sample", "pho_sub")], phog_kin.m[,c("sample", "phog_kin")], all = T)
        df1$technology <- "global"
        df1 <- merge(df1, clinical[, c("Specimen.Label", "tumor_normal")], by.x = c("sample"), by.y = c("Specimen.Label"), all.x = T)
        
        df2 <- merge(pp_cancer_pep.m[,c("Participant_ID", "value")], phog_kin.m[phog_kin.m$tumor_normal == "Tumor",c("Participant.ID", "phog_kin")], by.y = c("Participant.ID"), by.x = c("Participant_ID"), all = T)
        df2$value[df2$value == "#N/A"] <- NA
        df2$value <- as.numeric(as.vector(df2$value))
        df2$value[is.na(df2$value) & !is.na(df2$phosphopeptide.heavy.ratio)] <- as.numeric(df2$phosphopeptide.heavy.ratio[is.na(df2$value) & !is.na(df2$phosphopeptide.heavy.ratio)])
        colnames(df2) <- c("Participant_ID", "pho_sub", "phog_kin")
        df2$technology <- "SRM"
        df2$tumor_normal <- "Tumor"
        col_tmp <- intersect(colnames(df1), colnames(df2))
        df <- rbind(df1[,col_tmp], df2[,col_tmp])
        df <- df[!is.na(df$tumor_normal),]
        df$tumor_normal <- factor(df$tumor_normal, levels = c("Tumor", "Normal"))
        
        p = ggplot(df[!is.na(remove_outliers(df$phog_kin, out_thres = 1.5)),])
        p = p + geom_point(aes(x=phog_kin, y=pho_sub, color = tumor_normal), alpha=0.9, stroke = 0)
        p <- p + facet_grid(technology~., scales = "free_y", space = "fixed")
        p = p + labs(x = paste0("global phosphorylaion abundance for ", kinase),
                     y=paste0("SRM/global phosphorylation abundance for ", sub, " ", rsd))
        p = p + theme_nogrid()
        p = p + theme(axis.title = element_text(size=10),
                      axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5),
                      axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
        p = p + theme(title = element_text(size = 8))
        p
        fn = paste0(subdir, cancer, "_", kinase, "_", sub, "_", rsd ,"_global~SRM.pdf")
        ggsave(file=fn, height=6, width=5, useDingbats=FALSE)
      }
    }
  }
}


# cis --------------------------------------------------------------------
colnames(trans_tab)[5] <- "peptide"
common_columns <- intersect(colnames(trans_tab), colnames(cis_tab))
tab <- rbind(trans_tab[, common_columns], cis_tab[, common_columns])
#for (cancer in c("CO", "OV", "BRCA")) {
for (cancer in c("CO")) {
  pp_cancer <- pp_tab[, c("Peptide.Sequence", "Protein.Name", colnames(pp_tab)[grepl(pattern = substr(cancer, start = 1, stop = 2), x = colnames(pp_tab))])]
  tab_cancer <- tab[tab$cancer_type ==  cancer,]
  
  ## input global data
  pro <- fread(input = paste(inputD, cancer,"/",prefix[cancer], "_PRO_formatted_normalized_replicate_averaged.txt",sep=""), data.table = F)
  
  pho <- fread(input = paste(inputD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_replicate_averaged.txt",sep=""), data.table = F)
  pho <- cbind(formatPhosphosite(phosphosite_vector = pho$Phosphosite, pho$Gene), pho[,!(colnames(pho) %in% c("Gene", "Phosphosite"))])
  
  resultDnow <- makeOutDir()
  subdir <- paste0(resultDnow, cancer, "/")
  dir.create(subdir)
  
  subdir1 <- paste0(subdir, "cis", "/")
  dir.create(subdir1)
  for (peptide in unique(pp_cancer$Peptide.Sequence)) {
    pp_cancer_pep <- pp_cancer[pp_cancer$Peptide.Sequence == toupper(peptide),]
    tab_cancer_pep <- tab_cancer[grepl(pattern = peptide, x = tab_cancer$peptide, ignore.case = T),]
    pp_cancer_pep.m <- melt(pp_cancer_pep, id.vars = c("Peptide.Sequence", "Protein.Name"))
    pp_cancer_pep.m$Participant_ID <- str_split_fixed(string = pp_cancer_pep.m$variable, pattern = "X", 2)[,2]
    
    sub <- unique(tab_cancer_pep$SUBSTRATE)
    rsd <- unique(tab_cancer_pep$SUB_MOD_RSD)
    pho_sub <- pho[pho$SUBSTRATE == sub & pho$SUB_MOD_RSD == rsd,]
    if (nrow(pho_sub) > 0 ){
      pho_sub.m <- melt(pho_sub, id.vars = c("SUBSTRATE", "transcript", "SUB_MOD_RSD"))
      colnames(pho_sub.m) <- c("substrate", "transcript", "SUB_MOD_RSD", "sample", "pho_sub")
      kinases <- unique(tab_cancer_pep$KINASE)
      kinases <- kinases[kinases!=sub]
      for (kinase in kinases) {
        phog_kin <- phog[phog$Gene == kinase,]
        phog_kin.m <- melt(phog_kin, id.vars = c("Gene"))
        colnames(phog_kin.m) <- c("kinase", "sample", "phog_kin")
        phog_kin.m <- merge(phog_kin.m, clinical[, c("Specimen.Label", "Participant.ID", "tumor_normal")], by.x = c("sample"), by.y = c("Specimen.Label"), all.x = T)
        df1 <- merge(pho_sub.m[,c("sample", "pho_sub")], phog_kin.m[,c("sample", "phog_kin")], all = T)
        df1$technology <- "global"
        df1 <- merge(df1, clinical[, c("Specimen.Label", "tumor_normal")], by.x = c("sample"), by.y = c("Specimen.Label"), all.x = T)
        
        df2 <- merge(pp_cancer_pep.m[,c("Participant_ID", "value")], phog_kin.m[phog_kin.m$tumor_normal == "Tumor",c("Participant.ID", "phog_kin")], by.y = c("Participant.ID"), by.x = c("Participant_ID"), all = T)
        df2$value[df2$value == "#N/A"] <- NA
        df2$value <- as.numeric(as.vector(df2$value))
        df2$value[is.na(df2$value) & !is.na(df2$phosphopeptide.heavy.ratio)] <- as.numeric(df2$phosphopeptide.heavy.ratio[is.na(df2$value) & !is.na(df2$phosphopeptide.heavy.ratio)])
        colnames(df2) <- c("Participant_ID", "pho_sub", "phog_kin")
        df2$technology <- "SRM"
        df2$tumor_normal <- "Tumor"
        col_tmp <- intersect(colnames(df1), colnames(df2))
        df <- rbind(df1[,col_tmp], df2[,col_tmp])
        df <- df[!is.na(df$tumor_normal),]
        df$tumor_normal <- factor(df$tumor_normal, levels = c("Tumor", "Normal"))
        
        p = ggplot(df[!is.na(remove_outliers(df$phog_kin, out_thres = 1.5)),])
        p = p + geom_point(aes(x=phog_kin, y=pho_sub, color = tumor_normal), alpha=0.9, stroke = 0)
        p <- p + facet_grid(technology~., scales = "free_y", space = "fixed")
        p = p + labs(x = paste0("global phosphorylaion abundance for ", kinase),
                     y=paste0("SRM/global phosphorylation abundance for ", sub, " ", rsd))
        p = p + theme_nogrid()
        p = p + theme(axis.title = element_text(size=10),
                      axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5),
                      axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
        p = p + theme(title = element_text(size = 8))
        p
        fn = paste0(subdir1, cancer, "_", kinase, "_", sub, "_", rsd ,"_global~SRM.pdf")
        ggsave(file=fn, height=6, width=5, useDingbats=FALSE)
      }
    }
  }
}

library(GGally)
ggpairs(df, aes(colour = tumor_normal, alpha = 0.4), columns = 1:2)





