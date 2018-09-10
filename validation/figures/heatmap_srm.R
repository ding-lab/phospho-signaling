# Yige Wu @ WashU 2018 Apr
## make heatmap showing the ks outlier samples really agree regarding the order of abundance


# source ------------------------------------------------------------------
library(readxl)
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source("https://bioconductor.org/biocLite.R")

# inputs ------------------------------------------------------------------
## input the date of receiving the data
email_date <- "2018-04-20"

## input the sample ID mapping info
clinical <- loadSampMap()

## input gene symbol to gene symbol in crosstab
protein_names <- unique(up_cancer$Protein.Name)
names(protein_names) <- c("RPS6KA3", "ERBB2", "AKT1", "MET", "NUP153", "PDGFRA", "MAPK3", "EGFR", "PRKDC", "MAPK1", "MAP2K1", "SRC", "GAB2", "MAZL1", "GTF2I", "GAB1")
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

# cis just SRM samples --------------------------------------------------------------------
colnames(trans_tab)[5] <- "peptide"
common_columns <- intersect(colnames(trans_tab), colnames(cis_tab))
tab <- rbind(trans_tab[, common_columns], cis_tab[, common_columns])
for (cancer in "OV") {
  ## input cross tab for SRM values
  pp_cancer <- pp_tab[, c("Peptide.Sequence", "Protein.Name", colnames(pp_tab)[grepl(pattern = substr(cancer, start = 1, stop = 2), x = colnames(pp_tab))])]
  up_cancer <- up_tab[, c("Peptide.Sequence", "Protein.Name", colnames(pp_tab)[grepl(pattern = substr(cancer, start = 1, stop = 2), x = colnames(pp_tab))])]
  
  tab_cancer <- tab[tab$cancer_type ==  cancer,]
  tab_cancer$SELF <- ifelse(as.vector(tab_cancer$KINASE) == as.vector(tab_cancer$SUBSTRATE), "cis",  "trans")
  
  samples <- unique(tab_cancer[, c("Specimen_Label", "Participant_ID")])
  samples <- merge(samples, clinical[, c("Specimen.Label", "Participant.ID", "tumor_normal")], by.x = c("Specimen_Label"), by.y = c("Specimen.Label"), all.x = T)
  
  ## input global data
  phog <- loadPhosphoproteinNormalizedTumor(cancer = cancer)
  pho <- loadPhosphositeNormalizedTumor(cancer = cancer)
  pho <- cbind(formatPhosphosite(phosphosite_vector = pho$Phosphosite, pho$Gene), pho[,!(colnames(pho) %in% c("Gene", "Phosphosite"))])
  pro <- loadProteinNormalizedTumor(cancer = cancer)
  
  subdir <- paste0(makeOutDir(resultD = resultD), cancer, "_justSRM_samples/")
  dir.create(subdir)
  
  for (peptide in unique(pp_cancer$Peptide.Sequence)[4]) {
    pp_cancer_pep <- pp_cancer[pp_cancer$Peptide.Sequence == toupper(peptide),]
    pp_cancer_pep.m <- melt(pp_cancer_pep, id.vars = c("Peptide.Sequence", "Protein.Name"))
    pp_cancer_pep.m$Participant_ID <- str_split_fixed(string = pp_cancer_pep.m$variable, pattern = "X", 2)[,2]

    for (self in c("cis")) {
      tab_pep <- tab[grepl(pattern = peptide, x = tab$peptide, ignore.case = T),]
      tab_cancer_pep <- tab_cancer[grepl(pattern = peptide, x = tab_cancer$peptide, ignore.case = T),]
      
      if (nrow(tab_pep) > 0) {
        
        sub <- unique(tab_pep$SUBSTRATE)
        rsd <- unique(tab_pep$SUB_MOD_RSD)
        pho_sub <- pho[pho$SUBSTRATE == sub & pho$SUB_MOD_RSD == rsd,]

        if (nrow(pho_sub) > 0){
          pho_sub.m <- melt(pho_sub, id.vars = c("SUBSTRATE", "transcript", "SUB_MOD_RSD"))
          colnames(pho_sub.m) <- c("substrate", "transcript", "SUB_MOD_RSD", "sample", "pho_sub")

          for (kinase in sub) {
            pro_kin <- pro[pro$Gene == kinase,]
            pro_kin.m <- melt(pro_kin, id.vars = c("Gene"))
            colnames(pro_kin.m) <- c("kinase", "sample", "pro_kin")
            
            df <- merge(pho_sub.m, pro_kin.m, all = T)
            df$Participant_ID <- sampID2partID(sampleID_vector = as.vector(df$sample), sample_map = clinical)
            
            df <- merge(df, pp_cancer_pep.m, by = c("Participant_ID"), all = T)
            df$value[df$value == "#N/A"] <- NA
            df$value <- as.numeric(as.vector(df$value))
            colnames(df)[colnames(df) == "value"] <- c("pho_sub.srm")
            
            df$kinase <- kinase
            df <- merge(df, tab_cancer_pep, by = c("Participant_ID"), all.x = T)
            
            df <- df[!is.na(df$pho_sub.srm), ]
            
            if (nrow(df) > 0 & length(unique(df$sample)) > 1) {
              df0 <- df
              for (upeptide in up_cancer$Peptide.Sequence[up_cancer$Protein.Name == pp_cancer_pep$Protein.Name]) {
                up_cancer_pep_sub <- up_cancer[up_cancer$Peptide.Sequence == upeptide,]
                up_cancer_pep_sub.m <- melt(up_cancer_pep_sub, id.vars = c("Peptide.Sequence", "Protein.Name"))
                up_cancer_pep_sub.m$Participant_ID <- str_split_fixed(string = up_cancer_pep_sub.m$variable, pattern = "X", 2)[,2]

                df <- merge(df0, up_cancer_pep_sub.m[,c("Participant_ID", "value")], by.x = c("Participant_ID"), by.y = c("Participant_ID"), all = T)
                df$value[df$value == "#N/A"] <- NA
                df$value <- as.numeric(as.vector(df$value))
                colnames(df)[colnames(df) == "value"] <- c("pro_kin.srm")
                
                df <- df[!is.na(df$sample),]
                df$sample <- as.vector(df$sample)
                df <- df[!duplicated(df$pro_kin),]
                
                cap <- min(max(abs(df$pro_kin), na.rm = T), 2)
                df$pro_kin_capped <- df$pro_kin
                df$pro_kin_capped[df$pro_kin > cap] <- cap
                df$pro_kin_capped[df$pro_kin < (-cap)] <- (-cap)
                df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$pro_kin)])
                
                p1 = ggplot(df)
                p1 = p1 + geom_tile(aes(x=sample, y= paste0(kinase, " protein\n(global)"), fill= pro_kin_capped, width=0.7, height=0.7), size=0.5)
                p1 = p1 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
                p1 = p1 + scale_fill_gradientn(name= "global kinase phosphorylation abundance", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
                p1 = p1 + theme_bw() + theme_nogrid()
                p1 = p1 + theme(axis.text.x = element_blank())
                p1 = p1 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5, face = "bold"))
                p1 = p1 + theme(axis.ticks=element_blank(), 
                                legend.text = element_text(size = 2), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
                p1 = p1 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
                p1 = p1 + guides(color = FALSE)
                p1 = p1 + theme(legend.key.height = unit(0.15, "line"))
                p1
        
                cap <- min(max(abs(df$pho_sub), na.rm = T),2)
                df$pho_sub_capped <- df$pho_sub
                df$pho_sub_capped[df$pho_sub > cap] <- cap
                df$pho_sub_capped[df$pho_sub < (-cap)] <- (-cap)
                df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$pro_kin)])
                
                p2 = ggplot(df)
                p2 = p2 + geom_tile(aes(x=sample, y= paste0(sub, ":", rsd, "\n(global)"), fill= pho_sub_capped, width=0.7, height=0.7), size=0.5)
                p2 = p2 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
                p2 = p2 + scale_fill_gradientn(name= "global substrate phosphosite abundance", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
                p2 = p2 + theme_bw() + theme_nogrid()
                p2 = p2 + theme(axis.text.x = element_blank())
                p2 = p2 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5, face = "bold"))
                p2 = p2 + theme(axis.ticks=element_blank(), 
                                legend.text = element_text(size = 2), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
                p2 = p2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
                p2 = p2 + guides(color = FALSE)
                p2 = p2 + theme(legend.key.height = unit(0.15, "line"))
                
                gb1 <- ggplot_build(p1)
                g <- ggplot_gtable(gb1)
                n1 <- length(gb1$layout$panel_ranges[[1]]$y.labels)
                
                gb2 <- ggplot_build(p2)
                gB <- ggplot_gtable(gb2)
                g <- gtable:::rbind_gtable(g, gB, "last")
                n2 <- length(gb2$layout$panel_ranges[[1]]$y.labels)
                
                
                if(length(df$pro_kin.srm[!is.na(df$pro_kin.srm)]) > 0){
                  high <- quantile(x = df$pro_kin.srm, probs = 0.9, na.rm = T)
                  low <- quantile(x = df$pro_kin.srm, probs = 0.1, na.rm = T)
                  df$pro_kin.srm_capped <- df$pro_kin.srm
                  df$pro_kin.srm_capped[df$pro_kin.srm_capped > high] <- high
                  df$pro_kin.srm_capped[df$pro_kin.srm_capped < low] <- low
                  df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$pro_kin)])
                  
                  p3 = ggplot(df)
                  p3 = p3 + geom_tile(aes(x=sample, y= paste0(kinase, " protein\n(SRM)"), fill= pro_kin.srm_capped, width=0.7, height=0.7), size=0.5)
                  p3 = p3 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
                  p3 = p3 + scale_fill_gradientn(name= "SRM unmodified peptide abundance", na.value=NA, colours=RdBu1024, limit = c(low, high))
                  p3 = p3 + theme_bw() + theme_nogrid()
                  p3 = p3 + theme(axis.text.x = element_blank())
                  p3 = p3 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5, face = "bold"))
                  p3 = p3 + theme(axis.ticks=element_blank(), 
                                  legend.text = element_text(size = 2), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
                  p3 = p3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
                  p3 = p3 + guides(color = FALSE)
                  p3 = p3 + theme(legend.key.height = unit(0.15, "line"))
                  
                  gb3 <- ggplot_build(p3)
                  gC <- ggplot_gtable(gb3)
                  g <- gtable:::rbind_gtable(g, gC, "last")
                  n3 <- length(gb3$layout$panel_ranges[[1]]$y.labels)
                }

                if(length(df$pho_sub.srm[!is.na(df$pho_sub.srm)]) > 0){
                  high <- quantile(x = df$pho_sub.srm, probs = 0.9, na.rm = T)
                  low <- quantile(x = df$pho_sub.srm, probs = 0.1, na.rm = T)
                  df$pho_sub.srm_capped <- df$pho_sub.srm
                  df$pho_sub.srm_capped[df$pho_sub.srm_capped > high] <- high
                  df$pho_sub.srm_capped[df$pho_sub.srm_capped < low] <- low
                  df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$pro_kin)])
                  
                  p4 = ggplot(df)
                  p4 = p4 + geom_tile(aes(x=sample, y= paste0(sub, ":", rsd, "\n(SRM)"), fill= pho_sub.srm, width=0.7, height=0.7), size=0.5)
                  p4 = p4 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
                  p4 = p4 + scale_fill_gradientn(name= "SRM phosphopeptide abundance", na.value=NA, colours=RdBu1024, 
                                                 limit=c(min(df$pho_sub.srm, na.rm = T), 
                                                         max(df$pho_sub.srm, na.rm = T)))
                  p4 = p4 + theme_bw() + theme_nogrid()
                  p4 = p4 + theme(axis.text.x = element_blank())
                  p4 = p4 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5, face = "bold"))
                  p4 = p4 + theme(axis.ticks=element_blank(), 
                                  legend.text = element_text(size = 2), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
                  p4 = p4 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
                  p4 = p4 + guides(color = FALSE)
                  p4 = p4 + theme(legend.key.height = unit(0.15, "line"))
                  
                  gb4 <- ggplot_build(p4)
                  gD <- ggplot_gtable(gb4)
                  g <- gtable:::rbind_gtable(g, gD, "last")
                  n4 <- length(gb4$layout$panel_ranges[[1]]$y.labels)
                }
                panels <- g$layout$t[grep("panel", g$layout$name)]

                ## adjust height of each panel
                g$heights[panels] <- unit(x = rep(1, length(panels)), units = "null")

                subdir1 <- paste0(subdir, kinase, "_", sub, "_", rsd, "/")
                dir.create(subdir1)
                fn = paste0(subdir1, cancer, "_", kinase, "_", upeptide, "_", sub, "_", rsd, "_", peptide ,"_global~SRM.pdf")
                grid.newpage()
                pdf(fn, height = 1.5, width = 6)
                grid.draw(g)
                dev.off()
                
              }
            }
          }
        }
      }
    }
  }
}

# trans just SRM samples --------------------------------------------------------------------
colnames(trans_tab)[5] <- "peptide"
common_columns <- intersect(colnames(trans_tab), colnames(cis_tab))
tab <- rbind(trans_tab[, common_columns], cis_tab[, common_columns])
tab$SELF <- ifelse(as.vector(tab$KINASE) == as.vector(tab$SUBSTRATE), "cis",  "trans")
for (cancer in cancers_sort) {
  ## input cross tab for SRM values
  pp_cancer <- pp_tab[, c("Peptide.Sequence", "Protein.Name", colnames(pp_tab)[grepl(pattern = substr(cancer, start = 1, stop = 2), x = colnames(pp_tab))])]
  up_cancer <- up_tab[, c("Peptide.Sequence", "Protein.Name", colnames(pp_tab)[grepl(pattern = substr(cancer, start = 1, stop = 2), x = colnames(pp_tab))])]
  
  tab_cancer <- tab[tab$cancer_type ==  cancer,]
  
  samples <- unique(tab_cancer[, c("Specimen_Label", "Participant_ID")])
  samples <- merge(samples, clinical[, c("Specimen.Label", "Participant.ID", "tumor_normal")], by.x = c("Specimen_Label"), by.y = c("Specimen.Label"), all.x = T)
  
  ## input global data
  phog <- loadPhosphoproteinNormalizedTumor(cancer = cancer)
  pho <- loadPhosphositeNormalizedTumor(cancer = cancer)
  pho <- cbind(formatPhosphosite(phosphosite_vector = pho$Phosphosite, pho$Gene), pho[,!(colnames(pho) %in% c("Gene", "Phosphosite"))])
  pro <- loadProteinNormalizedTumor(cancer = cancer)
  
  subdir <- paste0(makeOutDir(resultD = resultD), cancer, "_justSRM_samples/")
  dir.create(subdir)
  
  for (peptide in unique(pp_cancer$Peptide.Sequence)) {
    pp_cancer_pep <- pp_cancer[pp_cancer$Peptide.Sequence == toupper(peptide),]
    pp_cancer_pep.m <- melt(pp_cancer_pep, id.vars = c("Peptide.Sequence", "Protein.Name"))
    pp_cancer_pep.m$Participant_ID <- str_split_fixed(string = pp_cancer_pep.m$variable, pattern = "X", 2)[,2]
    
    for (self in c("trans")) {
      tab_pep <- tab[grepl(pattern = peptide, x = tab$peptide, ignore.case = T) & tab$SELF == self,]
      tab_cancer_pep <- tab_cancer[grepl(pattern = peptide, x = tab_cancer$peptide, ignore.case = T),]
      
      if (nrow(tab_pep) > 0) {
        sub <- unique(tab_pep$SUBSTRATE)
        rsd <- unique(tab_pep$SUB_MOD_RSD)
        
        pho_sub <- pho[pho$SUBSTRATE == sub & pho$SUB_MOD_RSD == rsd,]
        if (nrow(pho_sub) > 0){
          pho_sub.m <- melt(pho_sub, id.vars = c("SUBSTRATE", "transcript", "SUB_MOD_RSD"))
          colnames(pho_sub.m) <- c("substrate", "transcript", "SUB_MOD_RSD", "sample", "pho_sub")
          
          kinases <- unique(tab_pep$KINASE)
          
          for (kinase in kinases) {
            pro_kin <- pro[pro$Gene == kinase,]
            pro_kin.m <- melt(pro_kin, id.vars = c("Gene"))
            colnames(pro_kin.m) <- c("kinase", "sample", "pro_kin")
            
            pro_sub <- pro[pro$Gene == sub,]
            pro_sub.m <- melt(pro_sub, id.vars = c("Gene"))
            colnames(pro_sub.m) <- c("sub", "sample", "pro_sub")
            
            phog_kin <- phog[phog$Gene == kinase,]
            phog_kin.m <- melt(phog_kin, id.vars = c("Gene"))
            colnames(phog_kin.m) <- c("kinase", "sample", "phog_kin")
            
            df <- merge(pho_sub.m, pro_kin.m, all = T)
            df$Participant_ID <- sampID2partID(sampleID_vector = as.vector(df$sample), sample_map = clinical)
            
            df <- merge(df, phog_kin.m, all = T)
            df <- merge(df, pro_sub.m, all = T)
            
            df <- merge(df, pp_cancer_pep.m, by = c("Participant_ID"), all = T)
            df$value[df$value == "#N/A"] <- NA
            df$value <- as.numeric(as.vector(df$value))
            colnames(df)[colnames(df) == "value"] <- c("pho_sub.srm")
            
            df <- merge(df, tab_cancer_pep, by = c("Participant_ID"), all.x = T)
            df <- df[!is.na(df$pho_sub.srm), ]
            
            if (nrow(df) > 0 & length(unique(df$sample)) > 1) {
              df0 <- df
              for (upeptide in up_cancer$Peptide.Sequence[up_cancer$Protein.Name == pp_cancer_pep$Protein.Name]) {
                up_cancer_pep_sub <- up_cancer[up_cancer$Peptide.Sequence == upeptide,]
                up_cancer_pep_sub.m <- melt(up_cancer_pep_sub, id.vars = c("Peptide.Sequence", "Protein.Name"))
                up_cancer_pep_sub.m$Participant_ID <- str_split_fixed(string = up_cancer_pep_sub.m$variable, pattern = "X", 2)[,2]
                
                df <- merge(df0, up_cancer_pep_sub.m[,c("Participant_ID", "value")], by.x = c("Participant_ID"), by.y = c("Participant_ID"), all = T)
                df$value[df$value == "#N/A"] <- NA
                df$value <- as.numeric(as.vector(df$value))
                colnames(df)[colnames(df) == "value"] <- c("pro_sub.srm")
                
                df <- merge(df0, up_cancer_pep_sub.m[,c("Participant_ID", "value")], by.x = c("Participant_ID"), by.y = c("Participant_ID"), all = T)
                df$value[df$value == "#N/A"] <- NA
                df$value <- as.numeric(as.vector(df$value))
                colnames(df)[colnames(df) == "value"] <- c("pro_sub.srm")
                
                df1 <- df
                for (upeptide_kin in up_cancer$Peptide.Sequence[up_cancer$Protein.Name == protein_names[kinase]]) {
                  up_cancer_pep_kin <- up_cancer[up_cancer$Peptide.Sequence == upeptide_kin,]
                  up_cancer_pep_kin.m <- melt(up_cancer_pep_kin, id.vars = c("Peptide.Sequence", "Protein.Name"))
                  up_cancer_pep_kin.m$Participant_ID <- str_split_fixed(string = up_cancer_pep_kin.m$variable, pattern = "X", 2)[,2]
                  
                  df <- merge(df1, up_cancer_pep_kin.m[,c("Participant_ID", "value")], by.x = c("Participant_ID"), by.y = c("Participant_ID"), all = T)
                  df$value[df$value == "#N/A"] <- NA
                  df$value <- as.numeric(as.vector(df$value))
                  colnames(df)[colnames(df) == "value"] <- c("pro_kin.srm")
                  
                  df <- df[!is.na(df$sample),]
                  df$sample <- as.vector(df$sample)
                  df <- df[!duplicated(df$pro_kin),]
                  
                  cap <- min(max(abs(df$phog_kin), na.rm = T), 2)
                  df$phog_kin_capped <- df$phog_kin
                  df$phog_kin_capped[df$phog_kin > cap] <- cap
                  df$phog_kin_capped[df$phog_kin < (-cap)] <- (-cap)
                  df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$phog_kin)])
                  
                  p1 = ggplot(df)
                  p1 = p1 + geom_tile(aes(x=sample, y= paste0(kinase, "\nphosphorylation\n(global)"), fill= phog_kin_capped, width=0.7, height=0.7), size=0.5)
                  p1 = p1 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
                  p1 = p1 + scale_fill_gradientn(name= "global kinase phosphorylation abundance", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
                  p1 = p1 + theme_bw() + theme_nogrid()
                  p1 = p1 + theme(axis.text.x = element_blank())
                  p1 = p1 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5, face = "bold"))
                  p1 = p1 + theme(axis.ticks=element_blank(), 
                                  legend.text = element_text(size = 2), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
                  p1 = p1 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
                  p1 = p1 + guides(color = FALSE)
                  p1 = p1 + theme(legend.key.height = unit(0.15, "line"))
                  gb1 <- ggplot_build(p1)
                  g <- ggplot_gtable(gb1)
                  n1 <- length(gb1$layout$panel_ranges[[1]]$y.labels)
                  h_count <- 1
                  
                  if(length(df$pro_kin[!is.na(df$pro_kin)]) > 0){
                    cap <- min(max(abs(df$pro_kin), na.rm = T), 2)
                    df$pro_kin_capped <- df$pro_kin
                    df$pro_kin_capped[df$pro_kin > cap] <- cap
                    df$pro_kin_capped[df$pro_kin < (-cap)] <- (-cap)
                    df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$phog_kin)])
                    
                    p1_ = ggplot(df)
                    p1_ = p1_ + geom_tile(aes(x=sample, y= paste0(kinase, " protein\n(global)"), fill= pro_kin_capped, width=0.7, height=0.7), size=0.5)
                    p1_ = p1_ + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
                    p1_ = p1_ + scale_fill_gradientn(name= "global kinase phosphorylation abundance", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
                    p1_ = p1_ + theme_bw() + theme_nogrid()
                    p1_ = p1_ + theme(axis.text.x = element_blank())
                    p1_ = p1_ + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5, face = "bold"))
                    p1_ = p1_ + theme(axis.ticks=element_blank(), 
                                      legend.text = element_text(size = 2), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
                    p1_ = p1_ + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
                    p1_ = p1_ + guides(color = FALSE)
                    p1_ = p1_ + theme(legend.key.height = unit(0.15, "line"))
                    
                    gb1_ <- ggplot_build(p1_)
                    gA_ <- ggplot_gtable(gb1_)
                    g <- gtable:::rbind_gtable(g, gA_, "last")
                    n1_ <- length(gb1_$layout$panel_ranges[[1]]$y.labels)
                    h_count <- h_count + 1
                  }
                  
                  if(length(df$pro_sub[!is.na(df$pro_sub)]) > 0){
                    cap <- min(max(abs(df$pro_sub), na.rm = T), 2)
                    df$pro_sub_capped <- df$pro_sub
                    df$pro_sub_capped[df$pro_sub > cap] <- cap
                    df$pro_sub_capped[df$pro_sub < (-cap)] <- (-cap)
                    df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$phog_kin)])
                    
                    p2_ = ggplot(df)
                    p2_ = p2_ + geom_tile(aes(x=sample, y= paste0(sub, " protein\n(global)"), fill= pro_sub_capped, width=0.7, height=0.7), size=0.5)
                    p2_ = p2_ + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
                    p2_ = p2_ + scale_fill_gradientn(name= "global substrate protein abundance", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
                    p2_ = p2_ + theme_bw() + theme_nogrid()
                    p2_ = p2_ + theme(axis.text.x = element_blank())
                    p2_ = p2_ + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5, face = "bold"))
                    p2_ = p2_ + theme(axis.ticks=element_blank(), 
                                      legend.text = element_text(size = 2), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
                    p2_ = p2_ + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
                    p2_ = p2_ + guides(color = FALSE)
                    p2_ = p2_ + theme(legend.key.height = unit(0.15, "line"))
                    
                    gb1_ <- ggplot_build(p2_)
                    gB_ <- ggplot_gtable(gb1_)
                    g <- gtable:::rbind_gtable(g, gB_, "last")
                    h_count <- h_count + 1
                  }
                  
                  cap <- min(max(abs(df$pho_sub), na.rm = T),2)
                  df$pho_sub_capped <- df$pho_sub
                  df$pho_sub_capped[df$pho_sub > cap] <- cap
                  df$pho_sub_capped[df$pho_sub < (-cap)] <- (-cap)
                  df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$phog_kin)])
                  
                  p2 = ggplot(df)
                  p2 = p2 + geom_tile(aes(x=sample, y= paste0(sub, ":", rsd, "\n(global)"), fill= pho_sub_capped, width=0.7, height=0.7), size=0.5)
                  p2 = p2 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
                  p2 = p2 + scale_fill_gradientn(name= "global substrate phosphosite abundance", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
                  p2 = p2 + theme_bw() + theme_nogrid()
                  p2 = p2 + theme(axis.text.x = element_blank())
                  p2 = p2 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5, face = "bold"))
                  p2 = p2 + theme(axis.ticks=element_blank(), 
                                  legend.text = element_text(size = 2), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
                  p2 = p2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
                  p2 = p2 + guides(color = FALSE)
                  p2 = p2 + theme(legend.key.height = unit(0.15, "line"))
                  h_count <- h_count + 1
                  
                  gb2 <- ggplot_build(p2)
                  gB <- ggplot_gtable(gb2)
                  g <- gtable:::rbind_gtable(g, gB, "last")
                  n2 <- length(gb2$layout$panel_ranges[[1]]$y.labels)
                  
                  if(length(df$pro_kin.srm[!is.na(df$pro_kin.srm)]) > 0){
                    high <- quantile(x = df$pro_kin.srm, probs = 0.9, na.rm = T)
                    low <- quantile(x = df$pro_kin.srm, probs = 0.1, na.rm = T)
                    df$pro_kin.srm_capped <- df$pro_kin.srm
                    df$pro_kin.srm_capped[df$pro_kin.srm_capped > high] <- high
                    df$pro_kin.srm_capped[df$pro_kin.srm_capped < low] <- low
                    df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$phog_kin)])
                    
                    p3 = ggplot(df)
                    p3 = p3 + geom_tile(aes(x=sample, y= paste0(kinase, " protein\n(SRM)"), fill= pro_kin.srm_capped, width=0.7, height=0.7), size=0.5)
                    p3 = p3 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
                    p3 = p3 + scale_fill_gradientn(name= "SRM unmodified peptide abundance", na.value=NA, colours=RdBu1024, limit = c(low, high))
                    p3 = p3 + theme_bw() + theme_nogrid()
                    p3 = p3 + theme(axis.text.x = element_blank())
                    p3 = p3 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5, face = "bold"))
                    p3 = p3 + theme(axis.ticks=element_blank(), 
                                    legend.text = element_text(size = 2), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
                    p3 = p3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
                    p3 = p3 + guides(color = FALSE)
                    p3 = p3 + theme(legend.key.height = unit(0.15, "line"))
                    
                    gb3 <- ggplot_build(p3)
                    gC <- ggplot_gtable(gb3)
                    g <- gtable:::rbind_gtable(g, gC, "last")
                    n3 <- length(gb3$layout$panel_ranges[[1]]$y.labels)
                    h_count <- h_count + 1
                    
                  }
                  
                  if(length(df$pro_sub.srm[!is.na(df$pro_sub.srm)]) > 0){
                    high <- quantile(x = df$pro_sub.srm, probs = 0.9, na.rm = T)
                    low <- quantile(x = df$pro_sub.srm, probs = 0.1, na.rm = T)
                    df$pro_sub.srm_capped <- df$pro_sub.srm
                    df$pro_sub.srm_capped[df$pro_sub.srm_capped > high] <- high
                    df$pro_sub.srm_capped[df$pro_sub.srm_capped < low] <- low
                    df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$phog_kin)])
                    
                    p3_ = ggplot(df)
                    p3_ = p3_ + geom_tile(aes(x=sample, y= paste0(sub, " protein\n(SRM)"), fill= pro_sub.srm_capped, width=0.7, height=0.7), size=0.5)
                    p3_ = p3_ + scale_color_manual(values = c("up_regulated" = "#E3_1A1C", "down_regulated" = "#1F78B4"))
                    p3_ = p3_ + scale_fill_gradientn(name= "SRM unmodified peptide abundance", na.value=NA, colours=RdBu1024, limit = c(low, high))
                    p3_ = p3_ + theme_bw() + theme_nogrid()
                    p3_ = p3_ + theme(axis.text.x = element_blank())
                    p3_ = p3_ + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5, face = "bold"))
                    p3_ = p3_ + theme(axis.ticks=element_blank(), 
                                      legend.text = element_text(size = 2), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
                    p3_ = p3_ + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
                    p3_ = p3_ + guides(color = FALSE)
                    p3_ = p3_ + theme(legend.key.height = unit(0.15, "line"))
                    
                    gb3_ <- ggplot_build(p3_)
                    gC_ <- ggplot_gtable(gb3_)
                    g <- gtable:::rbind_gtable(g, gC_, "last")
                    n3_ <- length(gb3_$layout$panel_ranges[[1]]$y.labels)
                    h_count <- h_count + 1
                  }
                  
                  if(length(df$pho_sub.srm[!is.na(df$pho_sub.srm)]) > 0){
                    high <- quantile(x = df$pho_sub.srm, probs = 0.9, na.rm = T)
                    low <- quantile(x = df$pho_sub.srm, probs = 0.1, na.rm = T)
                    df$pho_sub.srm_capped <- df$pho_sub.srm
                    df$pho_sub.srm_capped[df$pho_sub.srm_capped > high] <- high
                    df$pho_sub.srm_capped[df$pho_sub.srm_capped < low] <- low
                    df$sample <- factor(df$sample, levels = as.vector(df$sample)[order(df$phog_kin)])
                    
                    p4 = ggplot(df)
                    p4 = p4 + geom_tile(aes(x=sample, y= paste0(sub, ":", rsd, "\n(SRM)"), fill= pho_sub.srm, width=0.7, height=0.7), size=0.5)
                    p4 = p4 + scale_color_manual(values = c("up_regulated" = "#E31A1C", "down_regulated" = "#1F78B4"))
                    p4 = p4 + scale_fill_gradientn(name= "SRM phosphopeptide abundance", na.value=NA, colours=RdBu1024, 
                                                   limit=c(min(df$pho_sub.srm, na.rm = T), 
                                                           max(df$pho_sub.srm, na.rm = T)))
                    p4 = p4 + theme_bw() + theme_nogrid()
                    p4 = p4 + theme(axis.text.x = element_blank())
                    p4 = p4 + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 5, face = "bold"))
                    p4 = p4 + theme(axis.ticks=element_blank(), 
                                    legend.text = element_text(size = 2), legend.key.size = unit(1, "lines"), legend.title = element_text(size = 5))
                    p4 = p4 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
                    p4 = p4 + guides(color = FALSE)
                    p4 = p4 + theme(legend.key.height = unit(0.15, "line"))
                    
                    gb4 <- ggplot_build(p4)
                    gD <- ggplot_gtable(gb4)
                    g <- gtable:::rbind_gtable(g, gD, "last")
                    n4 <- length(gb4$layout$panel_ranges[[1]]$y.labels)
                    
                  }
                  panels <- g$layout$t[grep("panel", g$layout$name)]
                  
                  ## adjust height of each panel
                  g$heights[panels] <- unit(x = rep(1, length(panels)), units = "null")
                  
                  subdir1 <- paste0(subdir, kinase, "_", sub, "_", rsd, "/")
                  dir.create(subdir1)
                  fn = paste0(subdir1, cancer, "_", kinase, "_", upeptide_kin, "_", sub, "_", upeptide, "_", rsd, "_", peptide ,"_global~SRM.pdf")
                  grid.newpage()
                  pdf(fn, height = h_count*0.5, width = max(6*length(unique(df$sample))/15, 6))
                  grid.draw(g)
                  dev.off()
                  
                }
              }
            }
          }
        }
      }
    }
  }
}
