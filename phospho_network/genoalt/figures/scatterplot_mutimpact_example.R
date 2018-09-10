# Yige Wu @ WashU 2018 Apr
## check regulated pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R')

library(ggrepel)
library(readxl)
# set variables -----------------------------------------------------------
amp_log2cn <- 0.3
del_log2cn <- -0.3

# inputs ------------------------------------------------------------------
## input druggable gene list
clinical <- fread(input = paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180307.txt"), data.table = F)
subtype_map <- read_excel(paste0(cptac_sharedD, "5_CPTAC2_Breast_Prospective_Collection_BI/proteome-data-v1.01011/data/20171106_CPTAC2_ProspectiveBreastCancer_Sample_Annotations_Table_v79.xlsx"))
msi_score <- fread(input = "./cptac2p/analysis_results/phospho_network/genoalt/figures/scatterplot_mutimpact_example/CPTAC.MSI.score.tsv", data.table = F)
partIDs_msih <- msi_score$Sample[msi_score$Score >= 3.5]
partIDs_msil <- msi_score$Sample[msi_score$Score < 3.5 & msi_score$Score >= 1.0]
partIDs_mss <- msi_score$Sample[msi_score$Score < 1.0]

## input enzyme_substrate table
# ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD_extended.csv"))
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor/omnipath_networkin_DEPOD_SignorNotSiteMapped.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)

## input significant mutation impact
mut_cnv_sig_cans <- fread(input = paste0(ppnD, "genoalt/figures/boxplot_cnvimpact/mut_cnv_sig_cans.txt"), data.table = F)

## input kinases
all_kinases <- fread(input = "./Ding_Lab/Projects_Current/Gene Lists/All_Kinase_cleaned.txt", data.table = F)
# check every SMG of each cancer and overlap with kinase list ------------------------------------------
SMGs <- list()
CNVs <- list()

SMGs[["BRCA"]] <- c("PIK3CA", "TP53", "MAP3K1", "MAP2K4", "GATA3", "MLL3", "CDH1", "PTEN", "PIK3R1", "AKT1", "RUNX1", "CBFB", "TBX3", "NCOR1", "CTCF", "FOXA1", "SF3B1", "CDKN1B", "RB1", "AFF2", "NF1", "PTPN22", "PTPRD")
CNVs[["BRCA"]] <- c("PIK3CA", "ERBB2", "TP53", "MAP2K4", "MLL3", "CDKN2A", "PTEN", "RB1")

SMGs[["OV"]] <- c("TP53", "NF1", "BRCA1", "BRCA2", "RB1", "CDK12")

SMGs[["CO"]] <- c("APC", "TP53", "KRAS", "PIK3CA", "FBXW7", "SMAD4", "TCF7L2", "NRAS")

for (cancer in cancers2process) {
  ## see what's not covered in the consensus gene list
  print(SMGs[[cancer]][!(SMGs[[cancer]] %in% loadGeneList(gene_type = "driver", cancer = cancer, is.soft.limit = ""))])
  
  ## see how many SMGs per cancer types
  print(length(SMGs[[cancer]]))
  
  ## see how many SMGs are kinases
  print(length(SMGs[[cancer]][SMGs[[cancer]] %in% all_kinases$gene]))
  print(SMGs[[cancer]][SMGs[[cancer]] %in% all_kinases$gene])
}

# plot mutation impact -----------------------------------------------------------
resultDnow <- makeOutDir(resultD = resultD)
subdir1 <- paste0(resultDnow, "mutation", "/")
dir.create(subdir1)

for (cancer in c("BRCA")) {
  subdir2 <- paste0(subdir1, cancer, "/")
  dir.create(subdir2)
  
  ## input the complete table of mutation impact
  mut_cnv_tab <- fread(input = paste0(ppnD, "genoalt/tables/mut_cnv_impact/", cancer, "_mut_cnv_tab.txt"), data.table = F)
  mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE)
  mut_cnv_tab <- mut_cnv_tab[mut_cnv_tab$pair_pro %in% ptms_site_pairs_sup$pair_pro,]
  mut_cnv_tab$pair <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE, ":", mut_cnv_tab$SUB_MOD_RSD)
  mut_cnv_tab <- merge(mut_cnv_tab, unique(ptms_site_pairs_sup[, c("GENE", "enzyme_type")]), all.x = T)
  
  ## extract the SMG table
  # tab_mut <- mut_cnv_sig_cans[mut_cnv_sig_cans$cancer == cancer & mut_cnv_sig_cans$genoalt_type == "mutation",]
  # tab_mut <- tab_mut[(tab_mut$enzyme_driver_type == "oncogene" & tab_mut$meddiff > 0) | (tab_mut$enzyme_driver_type == "tsg" & tab_mut$meddiff < 0), ]
  tab_mut <- mut_cnv_tab[mut_cnv_tab$GENE %in% SMGs[[cancer]],]
  
  ## annoate substrate genes to TCGA pathways
  sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab_mut$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
  sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                   SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
  tab_mut <- merge(tab_mut, sub_genes2pathways, all.x = T)
  tab_mut <- tab_mut[!is.na(tab_mut$SUB_GENE.path),]
  
  pro_data <- loadProteinNormalizedTumor(cancer = cancer)
  pho_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
  phog_data <- loadPhosphoproteinNormalizedTumor(cancer = cancer)
  pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
  maf <- loadMaf(cancer = cancer, maf_files = maf_files)
  if (cancer %in% c("OV", "CO")) {
    cnv <- fread(input = paste0("./cptac2p/copy_number/gatk4wxscnv/v1.3.swap_contamination_fixed/prospective_somatic/deliverables/", tolower(cancer),"/gene_level_CNV.", tolower(cancer),".v1.3.2018-03-19.tsv"), data.table = F)
  }
  if (cancer %in% c("BRCA")) {
    cnv <- fread(input = "./cptac2p/copy_number/gatk4wxscnv/v1.3.swap_contamination_fixed/prospective_somatic/deliverables/br/gene_level_CNV.br.v1.3.2018-03-19.tsv", data.table = F)
  }
  
  for (i in 1:nrow(tab_mut)) {
    enzyme <- as.character(tab_mut[i, "GENE"])
    substrate <- as.character(tab_mut[i, "SUB_GENE"])
    rsd <- as.character(tab_mut[i, "SUB_MOD_RSD"])
    subdir3 <- paste0(subdir2, enzyme, "/")
    dir.create(subdir3)
    subdir4 <- paste0(subdir3, substrate, "/")
    dir.create(subdir4)
    subdir5 <- paste0(subdir4, rsd, "/")
    dir.create(subdir5)
    
    ## input phospho level
    pho_sub <- pho_data[pho_head$SUBSTRATE == substrate & pho_head$SUB_MOD_RSD == rsd,]

    pho.m <- melt(pho_sub)
    colnames(pho.m)[ncol(pho.m)] <- "pho_sub"

    sup_tab <- pho.m[, c("variable", "pho_sub")]
    sup_tab$partID <- sampID2partID(sampleID_vector = as.vector(sup_tab$variable), sample_map = clinical)
    
    ## input protein level of the enzyme
    pro_en <- pro_data[pro_data$Gene == enzyme,]
    pro_en.m <- melt(pro_en)
    colnames(pro_en.m) <- c("enzyme", "variable", "pro_en")
    sup_tab <- merge(sup_tab, pro_en.m[, c("variable", "pro_en")], all.x = T)
    
    ## input protein level of the substrate
    pro_sub <- pro_data[pro_data$Gene == substrate,]
    pro_sub.m <- melt(pro_sub)
    colnames(pro_sub.m) <- c("substrate", "variable", "pro_sub")
    sup_tab <- merge(sup_tab, pro_sub.m[, c("variable", "pro_sub")], all.x = T)
    
    ## input maf
    maf_en <- maf[(maf$Hugo_Symbol == enzyme),]
    # maf <- maf[(maf$Hugo_Symbol == enzyme | maf$Hugo_Symbol == substrate | maf$Hugo_Symbol %in% unique(ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$SUB_GENE == substrate & ptms_site_pairs_sup$SUB_MOD_RSD == rsd])),]
    maf_en <- maf_en[maf_en$Variant_Classification != "Silent",]
    if (nrow(maf_en) > 0) {
      maf_en$partID <- str_split_fixed(string = maf_en$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
      maf_en$aa_change <- paste0(maf_en$Hugo_Symbol, ":", maf_en$HGVSp_Short)
      maf_en$is.upstream <- ifelse(maf_en$Hugo_Symbol == enzyme, TRUE, FALSE)
      sup_tab <- merge(sup_tab, maf_en[, c("partID", "Variant_Classification", "aa_change", "is.upstream")], all.x = T) 
    }
    maf_sub <- maf[(maf$Hugo_Symbol == substrate),];  maf_sub <- maf_sub[maf_sub$Variant_Classification != "Silent",]
    if (nrow(maf_sub) > 0) {
      maf_sub$partID <- str_split_fixed(string = maf_sub$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
      maf_sub$aa_change <- paste0(maf_sub$Hugo_Symbol, ":", maf_sub$HGVSp_Short)
      sup_tab <- merge(sup_tab, maf_sub[, c("partID", "Variant_Classification", "aa_change")], by = c("partID"), all.x = T, suffixes = c("", ".sub")) 
    }
    
    cnv_en <- cnv[cnv$gene == enzyme,]
    cnv_en.m <- melt(cnv_en)
    colnames(cnv_en.m) <- c("enzyme", "partID", "log2cn")
    sup_tab <- merge(sup_tab, cnv_en.m[, c("partID", "log2cn")], all.x = T) 
    ## annotate sample type
    sup_tab$sample_type <- "control"
    sup_tab$sample_type[grepl(x = sup_tab$aa_change, pattern = enzyme)] <- "enzyme_mutation"
    sup_tab$sample_type[sup_tab$log2cn >= amp_log2cn] <- "enzyme_amplification"
    sup_tab$sample_type[sup_tab$log2cn <= del_log2cn] <- "enzyme_deletion"
    if (cancer == "BRCA") {
      sup_tab$subtype <- sampID2pam50(sampleID_vector = as.vector(sup_tab$variable), pam50_map = loadPAM50Map(), sample_map = loadSampMap())
    }
    if (cancer == "CO") {
      sup_tab$subtype <- "other"
      sup_tab$subtype[sup_tab$partID %in% partIDs_msih] <- "MSI-High"
      sup_tab$subtype[sup_tab$partID %in% partIDs_msil] <- "MSI-Low"
      sup_tab$subtype[sup_tab$partID %in% partIDs_mss] <- "MSS"
    }
    tab2p <- sup_tab
    tab2p <- unique(tab2p)
    tab2p$x <- tab2p$sample_type
    tab2p$x <- factor(tab2p$x, levels = c("enzyme_mutation", "enzyme_amplification", "enzyme_deletion", "control"))
    tab2p$text <- NA
    tab2p$text[tab2p$sample_type == "enzyme_mutation"] <- tab2p$aa_change[tab2p$sample_type == "enzyme_mutation"]
    tab2p$text[tab2p$sample_type == "enzyme_amplification" | tab2p$sample_type == "enzyme_deletion"] <- signif(tab2p$log2cn[tab2p$sample_type == "enzyme_amplification" | tab2p$sample_type == "enzyme_deletion"], digits = 2)
    tab2p$text[(tab2p$sample_type == "enzyme_amplification" | tab2p$sample_type == "enzyme_deletion") & !is.na(tab2p$aa_change)] <- paste0(tab2p$aa_change[(tab2p$sample_type == "enzyme_amplification" | tab2p$sample_type == "enzyme_deletion") & !is.na(tab2p$aa_change)], "|", signif(tab2p$log2cn[(tab2p$sample_type == "enzyme_amplification" | tab2p$sample_type == "enzyme_deletion") & !is.na(tab2p$aa_change)], digits = 2))
    tab2p$text[!is.na(tab2p$aa_change.sub) & !is.na(tab2p$text)] <- paste0(tab2p$text, "|", as.vector(tab2p$aa_change.sub[!is.na(tab2p$aa_change.sub) & is.na(tab2p$text)]))
    tab2p$text[!is.na(tab2p$aa_change.sub) & is.na(tab2p$text)] <- as.vector(tab2p$aa_change.sub[!is.na(tab2p$aa_change.sub) & is.na(tab2p$text)])
    
    pos <- position_jitter(width = 0.5, seed = 1)
    p = ggplot(tab2p, aes(x=x, y=pho_sub))
    p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
    p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
    p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type), 
                            force = 1, segment.size = 0.5, segment.alpha = 0.2, size=2.5,alpha=0.8, position = pos)
    p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_amplification", "enzyme_deletion", "control"),
                             label = c(paste0(enzyme, "_mutation"), paste0(enzyme, "_amplification"), paste0(enzyme, "_deletion"), "control"))
    p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio) in", cancer))
    p = p + theme_nogrid()
    p = p + theme(axis.title.x = element_blank(),
                  axis.text.x = element_text(face="bold", size=12, vjust=0.5, angle = 15),
                  axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
    p = p + theme(title = element_text(size = 8))
    p = p + scale_color_manual(values = c("enzyme_mutation" = "purple", "enzyme_amplification" = set1[1], "enzyme_deletion"  = set1[2], "control" = "grey"))
    p
    fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pho_sub.pdf")
    ggsave(file=fn, height=5, width=6, useDingbats=FALSE)
    
    pos <- position_jitter(width = 0.5, seed = 1)
    p = ggplot(tab2p, aes(x=x, y=pro_en))
    p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
    p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
    p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type), 
                            force = 1, segment.size = 0.5, segment.alpha = 0.2, size=2.5,alpha=0.8, position = pos)
    p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_amplification", "enzyme_deletion", "control"),
                             label = c(paste0(enzyme, "_mutation"), paste0(enzyme, "_amplification"), paste0(enzyme, "_deletion"), "control"))
    p = p + labs(y=paste0(enzyme, " protein abundance(log2 ratio) in", cancer))
    p = p + theme_nogrid()
    p = p + theme(axis.title.x = element_blank(),
                  axis.text.x = element_text(face="bold", size=12, vjust=0.5, angle = 15),
                  axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
    p = p + theme(title = element_text(size = 8))
    p = p + scale_color_manual(values = c("enzyme_mutation" = "purple", "enzyme_amplification" = set1[1], "enzyme_deletion"  = set1[2], "control" = "grey"))
    p
    fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_enzyme.pdf")
    ggsave(file=fn, height=5, width=6, useDingbats=FALSE)
    
    pos <- position_jitter(width = 0.5, seed = 1)
    p = ggplot(tab2p, aes(x=x, y=pro_sub))
    p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
    p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
    p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type), 
                            force = 1, segment.size = 0.5, segment.alpha = 0.2, size=2.5,alpha=0.8, position = pos)
    p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_amplification", "enzyme_deletion", "control"),
                             label = c(paste0(enzyme, "_mutation"), paste0(enzyme, "_amplification"), paste0(enzyme, "_deletion"), "control"))
    p = p + labs(y=paste0(substrate, " protein abundance(log2 ratio) in", cancer))
    p = p + theme_nogrid()
    p = p + theme(axis.title.x = element_blank(),
                  axis.text.x = element_text(face="bold", size=12, vjust=0.5, angle = 15),
                  axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
    p = p + theme(title = element_text(size = 8))
    p = p + scale_color_manual(values = c("enzyme_mutation" = "purple", "enzyme_amplification" = set1[1], "enzyme_deletion"  = set1[2], "control" = "grey"))
    p
    fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_sub.pdf")
    ggsave(file=fn, height=5, width=6, useDingbats=FALSE)
  }
}

# plot cnv impact -----------------------------------------------------------
resultDnow <- makeOutDir(resultD = resultD)

for (genoalt_type in c("amplification", "deletion")) {
  subdir1 <- paste0(resultDnow, genoalt_type, "/")
  dir.create(subdir1)
  for (cancer in cancers_sort) {
    tab_mut <- mut_cnv_sig_cans[mut_cnv_sig_cans$cancer == cancer & mut_cnv_sig_cans$genoalt_type == genoalt_type,]
    if (genoalt_type == "amplification") {
      tab_mut <- tab_mut[(tab_mut$enzyme_driver_type == "oncogene" & tab_mut$substrate_driver_type != "") | (tab_mut$GENE %in% ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "phosphatase"]), ]
    }
    if (genoalt_type == "deletion") {
      tab_mut <- tab_mut[(tab_mut$enzyme_driver_type == "tsg" & tab_mut$substrate_driver_type != "") | (tab_mut$GENE %in% ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "phosphatase"]), ]
    }
    pro_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PRO_formatted_normalized_noControl.txt",sep=""), data.table = F)
    pho_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
    phog_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt",sep=""), data.table = F)
    pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
    maf <- loadMaf(cancer = cancer, maf_files = maf_files)
    if (cancer %in% c("OV", "CO")) {
      cnv <- fread(input = paste0("./cptac2p/copy_number/gatk4wxscnv/v1.3.swap_contamination_fixed/prospective_somatic/deliverables/", tolower(cancer),"/gene_level_CNV.", tolower(cancer),".v1.3.2018-03-19.tsv"), data.table = F)
    }
    if (cancer %in% c("BRCA")) {
      cnv <- fread(input = "./cptac2p/copy_number/gatk4wxscnv/v1.3.swap_contamination_fixed/prospective_somatic/deliverables/br/gene_level_CNV.br.v1.3.2018-03-19.tsv", data.table = F)
    }
    
    for (i in 1:nrow(tab_mut)) {
      enzyme <- as.character(tab_mut[i, "GENE"])
      substrate <- as.character(tab_mut[i, "SUB_GENE"])
      rsd <- as.character(tab_mut[i, "SUB_MOD_RSD"])
      subdir2 <- paste0(subdir1, cancer, "/")
      dir.create(subdir2)
      subdir3 <- paste0(subdir2, enzyme, "/")
      dir.create(subdir3)
      subdir4 <- paste0(subdir3, substrate, "/")
      dir.create(subdir4)
      subdir5 <- paste0(subdir4, rsd, "/")
      dir.create(subdir5)
      
      ## input phospho level
      pho_sub <- pho_data[pho_head$SUBSTRATE == substrate & pho_head$SUB_MOD_RSD == rsd,]
      
      pho.m <- melt(pho_sub)
      colnames(pho.m)[ncol(pho.m)] <- "pho_sub"
      
      sup_tab <- pho.m[, c("variable", "pho_sub")]
      sup_tab$partID <- sampID2partID(sampleID_vector = as.vector(sup_tab$variable), sample_map = clinical)
      
      ## input protein level of the enzyme
      pro_en <- pro_data[pro_data$Gene == enzyme,]
      pro_en.m <- melt(pro_en)
      colnames(pro_en.m) <- c("enzyme", "variable", "pro_en")
      sup_tab <- merge(sup_tab, pro_en.m[, c("variable", "pro_en")], all.x = T)
      
      ## input protein level of the substrate
      pro_sub <- pro_data[pro_data$Gene == substrate,]
      pro_sub.m <- melt(pro_sub)
      colnames(pro_sub.m) <- c("substrate", "variable", "pro_sub")
      sup_tab <- merge(sup_tab, pro_sub.m[, c("variable", "pro_sub")], all.x = T)
      
      ## input maf
      maf_en <- maf[(maf$Hugo_Symbol == enzyme),];  maf_en <- maf_en[maf_en$Variant_Classification != "Silent",]
      if (nrow(maf_en) > 0) {
        maf_en$partID <- str_split_fixed(string = maf_en$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
        maf_en$aa_change <- paste0(maf_en$Hugo_Symbol, ":", maf_en$HGVSp_Short)
        maf_en$is.upstream <- ifelse(maf_en$Hugo_Symbol == enzyme, TRUE, FALSE)
        sup_tab <- merge(sup_tab, maf_en[, c("partID", "Variant_Classification", "aa_change", "is.upstream")], all.x = T) 
      }
      
      maf_sub <- maf[(maf$Hugo_Symbol == substrate),];  maf_sub <- maf_sub[maf_sub$Variant_Classification != "Silent",]
      if (nrow(maf_sub) > 0) {
        maf_sub$partID <- str_split_fixed(string = maf_sub$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
        maf_sub$aa_change <- paste0(maf_sub$Hugo_Symbol, ":", maf_sub$HGVSp_Short)
        sup_tab <- merge(sup_tab, maf_sub[, c("partID", "Variant_Classification", "aa_change")], by = c("partID"), all.x = T, suffixes = c("", ".sub")) 
      }
      cnv_en <- cnv[cnv$gene == enzyme,]
      cnv_en.m <- melt(cnv_en)
      colnames(cnv_en.m) <- c("enzyme", "partID", "log2cn")
      sup_tab <- merge(sup_tab, cnv_en.m[, c("partID", "log2cn")], all.x = T) 
      ## annotate sample type
      sup_tab$sample_type <- "control"
      sup_tab$sample_type[grepl(x = sup_tab$aa_change, pattern = enzyme)] <- "enzyme_mutation"
      sup_tab$sample_type[sup_tab$log2cn >= amp_log2cn] <- "enzyme_amplification"
      sup_tab$sample_type[sup_tab$log2cn <= del_log2cn] <- "enzyme_deletion"
      if (cancer == "BRCA") {
        sup_tab$subtype <- sampID2pam50(sampleID_vector = as.vector(sup_tab$variable), subtype_map = subtype_map)
      }
      if (cancer == "CO") {
        sup_tab$subtype <- "other"
        sup_tab$subtype[sup_tab$partID %in% partIDs_msih] <- "MSI-High"
        sup_tab$subtype[sup_tab$partID %in% partIDs_msil] <- "MSI-Low"
        sup_tab$subtype[sup_tab$partID %in% partIDs_mss] <- "MSS"
      }
      if (cancer == "OV") {
        sup_tab$subtype <- ""
      }
      tab2p <- sup_tab
      tab2p <- unique(tab2p)
      tab2p$x <- tab2p$sample_type
      tab2p$x <- factor(tab2p$x, levels = c("enzyme_mutation", "enzyme_amplification", "enzyme_deletion", "control"))
      tab2p$text[tab2p$sample_type == "enzyme_mutation"] <- tab2p$aa_change[tab2p$sample_type == "enzyme_mutation"]
      tab2p$text[tab2p$sample_type == "enzyme_amplification" | tab2p$sample_type == "enzyme_deletion"] <- signif(tab2p$log2cn[tab2p$sample_type == "enzyme_amplification" | tab2p$sample_type == "enzyme_deletion"], digits = 2)
      tab2p$text[(tab2p$sample_type == "enzyme_amplification" | tab2p$sample_type == "enzyme_deletion") & !is.na(tab2p$aa_change)] <- paste0(tab2p$aa_change[(tab2p$sample_type == "enzyme_amplification" | tab2p$sample_type == "enzyme_deletion") & !is.na(tab2p$aa_change)], "|", signif(tab2p$log2cn[(tab2p$sample_type == "enzyme_amplification" | tab2p$sample_type == "enzyme_deletion") & !is.na(tab2p$aa_change)], digits = 2))
      tab2p$text[!is.na(tab2p$aa_change.sub) & !is.na(tab2p$text)] <- paste0(tab2p$text, "|", as.vector(tab2p$aa_change.sub[!is.na(tab2p$aa_change.sub) & is.na(tab2p$text)]))
      tab2p$text[!is.na(tab2p$aa_change.sub) & is.na(tab2p$text)] <- as.vector(tab2p$aa_change.sub[!is.na(tab2p$aa_change.sub) & is.na(tab2p$text)])
      
      pos <- position_jitter(width = 0.5, seed = 1)
      p = ggplot(tab2p, aes(x=x, y=pho_sub))
      p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
      p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
      p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type), 
                              force = 1, segment.size = 0.5, segment.alpha = 0.2, size=2.5, alpha=0.8, position = pos)
      p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio) in ", cancer))
      p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_amplification", "enzyme_deletion", "control"),
                               label = c(paste0(enzyme, "_mutation"), paste0(enzyme, "_amplification"), paste0(enzyme, "_deletion"), "control"))
      p = p + theme_nogrid()
      p = p + theme(axis.title.x = element_blank(),
                    axis.text.x = element_text(face="bold", size=8, vjust=0.5, angle = 30),
                    axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
      p = p + theme(title = element_text(size = 8))
      p = p + scale_color_manual(values = c("enzyme_mutation" = "purple", "enzyme_amplification" = set1[1], "enzyme_deletion"  = set1[2], "control" = "grey"))
      p
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pho_sub.pdf")
      ggsave(file=fn, height=5, width=6, useDingbats=FALSE)
      
      p = ggplot(tab2p, aes(x=x, y=pro_en))
      p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
      p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
      p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type), 
                              force = 1, segment.size = 0.5, segment.alpha = 0.2, size=2.5,alpha=0.8, position = pos)
      p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_amplification", "enzyme_deletion", "control"),
                               label = c(paste0(enzyme, "_mutation"), paste0(enzyme, "_amplification"), paste0(enzyme, "_deletion"), "control"))
      p = p + labs(y=paste0(enzyme, " protein abundance(log2 ratio) in", cancer))
      p = p + theme_nogrid()
      p = p + theme(axis.title.x = element_blank(),
                    axis.text.x = element_text(face="bold", size=12, vjust=0.5, angle = 15),
                    axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
      p = p + theme(title = element_text(size = 8))
      p = p + scale_color_manual(values = c("enzyme_mutation" = "purple", "enzyme_amplification" = set1[1], "enzyme_deletion"  = set1[2], "control" = "grey"))
      p
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_enzyme.pdf")
      ggsave(file=fn, height=5, width=6, useDingbats=FALSE)
      
      p = ggplot(tab2p, aes(x=x, y=pro_sub))
      p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
      p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
      p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type), 
                              force = 1, segment.size = 0.5, segment.alpha = 0.2, size=2.5,alpha=0.8, position = pos)
      p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_amplification", "enzyme_deletion", "control"),
                               label = c(paste0(enzyme, "_mutation"), paste0(enzyme, "_amplification"), paste0(enzyme, "_deletion"), "control"))
      p = p + labs(y=paste0(substrate, " protein abundance(log2 ratio) in", cancer))
      p = p + theme_nogrid()
      p = p + theme(axis.title.x = element_blank(),
                    axis.text.x = element_text(face="bold", size=12, vjust=0.5, angle = 15),
                    axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
      p = p + theme(title = element_text(size = 8))
      p = p + scale_color_manual(values = c("enzyme_mutation" = "purple", "enzyme_amplification" = set1[1], "enzyme_deletion"  = set1[2], "control" = "grey"))
      p
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_sub.pdf")
      ggsave(file=fn, height=5, width=6, useDingbats=FALSE)
    }
  }
  
}


