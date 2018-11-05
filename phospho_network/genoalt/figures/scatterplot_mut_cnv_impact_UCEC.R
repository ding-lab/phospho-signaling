# Yige Wu @ WashU 2018 Apr
## check regulated pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/expression_matrices/expression_shared.R')

library(ggrepel)
library(readxl)
library(ggpubr)
# set variables -----------------------------------------------------------
sig_thres <- 0.05
networkin_thres <- 5
cancers2process <- cancers_sort
amp_thres <- 0.1
del_thres <- -0.1
num_genoalt_thres <- 4

my_comparisons <- list( enzyme_mutation = c("enzyme_mutation", "control"), 
                        enzyme_deep_amplification = c("enzyme_deep_amplification", "control"), 
                        enzyme_deep_deletion = c("enzyme_deep_deletion", "control"), 
                        enzyme_shallow_amplification = c("enzyme_shallow_amplification", "control"), 
                        enzyme_shallow_deletion = c("enzyme_shallow_deletion", "control"))


# inputs ------------------------------------------------------------------
## input druggable gene list
clinical <- fread(input = paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180307.txt"), data.table = F)
msi_score <- fread(input = "./cptac2p/analysis_results/phospho_network/genoalt/figures/scatterplot_mutimpact_example/CPTAC.MSI.score.tsv", data.table = F)
partIDs_msih <- msi_score$Sample[msi_score$Score >= 3.5]
partIDs_msil <- msi_score$Sample[msi_score$Score < 3.5 & msi_score$Score >= 1.0]
partIDs_mss <- msi_score$Sample[msi_score$Score < 1.0]

## input enzyme_substrate table
## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source == "NetKIN" | (ptms_site_pairs_sup$Source != "NetKIN" & ptms_site_pairs_sup$networkin_score >= networkin_thres),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)

## input kinases
all_kinases <- fread(input = "./Ding_Lab/Projects_Current/Gene Lists/All_Kinase_cleaned.txt", data.table = F)
# check every SMG of each cancer and overlap with kinase list ------------------------------------------

## input drive gene list
driver_genes <- read_excel("./Ding_Lab/Projects_Current/TCGA_data/gene_lists/mmc1.xlsx", 
                           sheet = "Table S1", skip = 3)
driver_genes <- data.frame(driver_genes)
oncogenes <- driver_genes$Gene[grepl(x = driver_genes$Tumor.suppressor.or.oncogene.prediction..by.20.20.., pattern = "oncogene")]
tsgs <- driver_genes$Gene[grepl(x = driver_genes$Tumor.suppressor.or.oncogene.prediction..by.20.20.., pattern = "tsg")]
oncogenes <- c(oncogenes, "CTNND1", "PIK3R1")
tsgs <- tsgs[!(tsgs %in% c("CDH1", "CTNND1", "PIK3R1"))]
oncogenes <- unique(oncogenes)
tsgs <- unique(tsgs)

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

for (cancer in "UCEC") {
  subdir2 <- paste0(subdir1, cancer, "/")
  dir.create(subdir2)
  
  ## input the complete table of mutation impact
  mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_cna_impact_UCEC/mut_cnv_sig_cans.txt"), data.table = F)
  mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE)
  mut_cnv_tab$pair <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE, ":", mut_cnv_tab$SUB_MOD_RSD)
  mut_cnv_tab <- merge(mut_cnv_tab, unique(ptms_site_pairs_sup[, c("GENE", "enzyme_type")]), all.x = T)
  
  ##filter out bad annotated enzyme type
  mut_cnv_tab <- mut_cnv_tab[!(mut_cnv_tab$GENE == "PTPRG" & mut_cnv_tab$enzyme_type == "kinase"),]
  
  ## annotate cis and trans
  mut_cnv_tab$SELF <- "trans"
  mut_cnv_tab$SELF[as.vector(mut_cnv_tab$GENE) == as.vector(mut_cnv_tab$SUB_GENE)] <- "cis"
  
  ## input somatic mutation matrix
  mut_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_mutation_matrix_UCEC/", cancer, "_somatic_mutation_in_enzyme_substrate.txt"), data.table = F)
  ## input CNA matrix
  deep_del_thresholded_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_cna_matrix/", cancer, "_deep_deletions_in_enzyme_substrate.txt"), data.table = F)
  deep_amp_thresholded_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_cna_matrix/", cancer, "_deep_amplifications_in_enzyme_substrate.txt"), data.table = F)
  shallow_del_thresholded_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_cna_matrix/", cancer, "_shallow_deletions_", del_thres, "_in_enzyme_substrate.txt"), data.table = F)
  shallow_amp_thresholded_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_cna_matrix/", cancer,"_shallow_amplifications_", amp_thres, "_in_enzyme_substrate.txt"), data.table = F)
  
  maf <- loadMaf(cancer = cancer, maf_files = maf_files)
  
  cna <- fread(input = paste0(cptac2pD, "copy_number/gistic/outputs/cptac2p.", tolower(substr(x = cancer, start = 1, stop = 2)), "/somatic/run_gistic2/0.25/0.99/armpeel1/somatic_clean/focal_data_by_genes.txt"),
               data.table = F)
  for (self in c("trans")) {
    # self <- "trans"
    ## extract the SMG table
    # tab_mut <- mut_cnv_tab[(mut_cnv_tab$GENE %in% SMGs[[cancer]] | mut_cnv_tab$cancer == "OV") & mut_cnv_tab$SELF == self,]
    tab_mut <- mut_cnv_tab[mut_cnv_tab$SUB_GENE == "TP53" & mut_cnv_tab$GENE == "TP53",]
    tab_mut <- tab_mut[!duplicated(tab_mut$SUB_MOD_RSD),]
    # tab_mut <- mut_cnv_tab[mut_cnv_tab$SUB_GENE == "AKT1" & mut_cnv_tab$GENE == "AKT1",]
    # tab_mut <- tab_mut[tab_mut$GENE %in% driver_genes$Gene,]
    # tab_mut <- mut_cnv_tab[mut_cnv_tab$SUB_GENE == "PIK3R1" & mut_cnv_tab$p_mut < sig_thres & mut_cnv_tab$p_mut > 0 & mut_cnv_tab$SELF == self,]
    # tab_mut <- mut_cnv_tab[mut_cnv_tab$SUB_GENE == "AKT1" & mut_cnv_tab$p_mut < sig_thres & mut_cnv_tab$p_mut > 0 & mut_cnv_tab$SELF == self,]
    
    ## annoate substrate genes to TCGA pathways
    sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab_mut$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
    sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                     SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
    tab_mut <- merge(tab_mut, sub_genes2pathways, all.x = T)
    # tab_mut <- tab_mut[!is.na(tab_mut$SUB_GENE.path),]
    
    pro_data <- loadProteinNormalizedTumor(cancer = cancer)
    pho_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
    phog_data <- loadPhosphoproteinNormalizedTumor(cancer = cancer)
    pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
    
    if (nrow(tab_mut) == 0) {
      print(paste0(cancer, "no data!"))
    } else {
      for (i in 1:nrow(tab_mut)) {
        enzyme <- as.character(tab_mut[i, "GENE"])
        # enzyme <- "AKT1"
        substrate <- as.character(tab_mut[i, "SUB_GENE"])
        # substrate <- "AKT1"
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
        
        ## input phosphoprotein level of the enzyme
        phog_en <- phog_data[phog_data$Gene == enzyme,]
        phog_en.m <- melt(phog_en)
        colnames(phog_en.m) <- c("enzyme", "variable", "phog_en")
        sup_tab <- merge(sup_tab, phog_en.m[, c("variable", "phog_en")], all.x = T)
        
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
        if (self == "trans") {
          maf_sub <- maf[(maf$Hugo_Symbol == substrate),];  maf_sub <- maf_sub[maf_sub$Variant_Classification != "Silent",]
          if (nrow(maf_sub) > 0) {
            maf_sub$partID <- str_split_fixed(string = maf_sub$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
            maf_sub$aa_change <- paste0(maf_sub$Hugo_Symbol, ":", maf_sub$HGVSp_Short)
            sup_tab <- merge(sup_tab, maf_sub[, c("partID", "Variant_Classification", "aa_change")], by = c("partID"), all.x = T, suffixes = c("", ".sub")) 
          }
        }
        
        cna_en <- cna[cna$`Gene Symbol` == enzyme,]
        cna_en.m <- melt(cna_en, id.vars = c("Gene Symbol", "Gene ID", "Cytoband"))
        colnames(cna_en.m) <- c("enzyme", "Gene ID", "Cytoband" , "partID", "log2cn")
        sup_tab <- merge(sup_tab, cna_en.m[, c("partID", "log2cn")], all.x = T)
        
        
        ## annotate sample type depending on the genomic alterations
        partIDs_overlap <- intersect(colnames(deep_del_thresholded_mat), intersect(unique(sup_tab$partID), colnames(mut_mat)))
        mut_mat_en <- mut_mat[mut_mat$Hugo_Symbol == enzyme,]
        if (nrow(mut_mat_en) > 0){
          mut_partIDs <- partIDs_overlap[(mut_mat_en[,partIDs_overlap] != "") & (mut_mat_en[, partIDs_overlap] != "Silent")]
        } else {
          mut_partIDs <- NULL 
        }
        deep_del_thresholded_mat_en <- deep_del_thresholded_mat[deep_del_thresholded_mat$Hugo_Symbol == enzyme,]
        if (nrow(deep_del_thresholded_mat_en) > 0){
          deep_del_partIDs <- partIDs_overlap[(deep_del_thresholded_mat_en[,partIDs_overlap] == "TRUE")]
        } else {
          deep_del_partIDs <- NULL 
        }
        deep_amp_thresholded_mat_en <- deep_amp_thresholded_mat[deep_amp_thresholded_mat$Hugo_Symbol == enzyme,]
        if (nrow(deep_amp_thresholded_mat_en) > 0){
          deep_amp_partIDs <- partIDs_overlap[(deep_amp_thresholded_mat_en[,partIDs_overlap] == "TRUE")]
        } else {
          deep_amp_partIDs <- NULL 
        }
        shallow_amp_thresholded_mat_en <- shallow_amp_thresholded_mat[shallow_amp_thresholded_mat$Hugo_Symbol == enzyme,]
        if (nrow(shallow_amp_thresholded_mat_en) > 0){
          shallow_amp_partIDs <- partIDs_overlap[(shallow_amp_thresholded_mat_en[,partIDs_overlap] == "TRUE")]
        } else {
          shallow_amp_partIDs <- NULL 
        }
        shallow_del_thresholded_mat_en <- shallow_del_thresholded_mat[shallow_del_thresholded_mat$Hugo_Symbol == enzyme,]
        if (nrow(shallow_del_thresholded_mat_en) > 0){
          shallow_del_partIDs <- partIDs_overlap[(shallow_del_thresholded_mat_en[,partIDs_overlap] == "TRUE")]
        } else {
          shallow_del_partIDs <- NULL 
        }
        
        sup_tab$sample_type <- "control"
        sup_tab$sample_type[sup_tab$partID %in% mut_partIDs] <- "enzyme_mutation"
        sup_tab$sample_type[sup_tab$partID %in% deep_del_partIDs] <- "enzyme_deep_deletion"
        sup_tab$sample_type[sup_tab$partID %in% deep_amp_partIDs] <- "enzyme_deep_amplification"
        sup_tab$sample_type[sup_tab$partID %in% shallow_del_partIDs] <- "enzyme_shallow_deletion"
        sup_tab$sample_type[sup_tab$partID %in% shallow_amp_partIDs] <- "enzyme_shallow_amplification"
        sample_type2test <- data.frame(table(unique(sup_tab[, c("sample_type", "partID")])[, "sample_type"]))
        sample_type2test <- as.vector(sample_type2test$Var1[sample_type2test$Freq >= num_genoalt_thres & sample_type2test$Var1 != "control"])
        
        sup_tab$subtype <- ""
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
        tab2p$x <- factor(tab2p$x, levels = c("enzyme_mutation", "enzyme_deep_amplification", "enzyme_deep_deletion", "enzyme_shallow_amplification", "enzyme_shallow_deletion", "control"))
        tab2p$text <- NA
        tab2p$text[tab2p$sample_type == "enzyme_mutation"] <- tab2p$aa_change[tab2p$sample_type == "enzyme_mutation"]
        tab2p$text[tab2p$sample_type != "enzyme_mutation" & tab2p$sample_type != "control"] <- signif(tab2p$log2cn[tab2p$sample_type != "enzyme_mutation" & tab2p$sample_type != "control"], digits = 2)
        row2edit <- (tab2p$sample_type != "enzyme_mutation") & !is.na(tab2p$aa_change)
        tab2p$text[row2edit] <- paste0(tab2p$aa_change[row2edit], "|", signif(tab2p$log2cn[row2edit]))
        row2edit <- (tab2p$sample_type == "enzyme_mutation") & (tab2p$partID %in% c(shallow_amp_partIDs, shallow_del_partIDs, deep_amp_partIDs, deep_del_partIDs))
        tab2p$text[row2edit] <- paste0(tab2p$aa_change[row2edit], "|", signif(tab2p$log2cn[row2edit]))
        
        
        if (self == "trans") {
          tab2p$text[!is.na(tab2p$aa_change.sub) & is.na(tab2p$text)] <- as.vector(tab2p$aa_change.sub[!is.na(tab2p$aa_change.sub) & is.na(tab2p$text)])
          tab2p$text[!is.na(tab2p$aa_change.sub) & !is.na(tab2p$text)] <- paste0(tab2p$text[!is.na(tab2p$aa_change.sub) & !is.na(tab2p$text)], "|",  tab2p$aa_change.sub[!is.na(tab2p$aa_change.sub) & !is.na(tab2p$text)])
        }
        
        
        # write table for protein paint -------------------------------------------
        if (self == "cis") {
          tab4jude_pho <- data.frame(text4proteinpaint = rsd, Position = strsplit(x = rsd, split = '[STY]')[[1]][2], judeClass = "M")
          tab4jude_mut <- data.frame(text4proteinpaint = maf_en$HGVSp_Short, 
                                     Position = str_split_fixed(string = maf_en$HGVSp_Short, pattern = '[A-Z]|[a-z]|\\*', n = 4)[, 3], 
                                     judeClass = substr(x = maf_en$Variant_Classification, start = 1, stop = 1))
          
          tab4jude <- rbind(tab4jude_pho, tab4jude_mut)
          write.table(x = tab4jude, file = paste0(subdir5, "tab4jude_", enzyme, "_", rsd, "_", cancer, ".txt"), row.names = F, quote = F, sep = ";", col.names = F)
        } else {
          tab4jude_mut <- data.frame(text4proteinpaint = maf_en$HGVSp_Short, 
                                     Position = str_split_fixed(string = maf_en$HGVSp_Short, pattern = '[A-Z]|[a-z]|\\*', n = 4)[, 3], 
                                     judeClass = substr(x = maf_en$Variant_Classification, start = 1, stop = 1))
          
          tab4jude <- tab4jude_mut
          write.table(x = tab4jude, file = paste0(subdir5, "tab4jude_", enzyme, "_", cancer, ".txt"), row.names = F, quote = F, sep = ";", col.names = F)
        }

        tab2p$y <- as.vector(tab2p$pho_sub)
        pos <- position_jitter(width = 0.5, seed = 1)
        p = ggplot(tab2p, aes(x=phog_en, y=y))
        p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
        p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                force = 1, segment.size = 0.5, segment.alpha = 0.2, size=2.5,alpha=0.8, position = pos)
        p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_deep_amplification", "enzyme_deep_deletion", "enzyme_shallow_amplification", "enzyme_shallow_deletion", "control"),
                                 label = c(paste0(enzyme, "_mutation"),
                                           paste0(enzyme, "_deep_amplification"), paste0(enzyme, "_deep_deletion"),
                                           paste0(enzyme, "_shallow_amplification"), paste0(enzyme, "_shallow_deletion"),
                                           "control"))
        p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio) in", cancer))
        p = p + theme_nogrid()
        p = p + theme(axis.title.x = element_blank(),
                      axis.text.x = element_text(size= 10, vjust=0.5, hjust = 0.5, angle = 90),
                      axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
        p = p + theme(title = element_text(size = 8))
        p = p + scale_color_manual(values = c("enzyme_mutation" = "purple",
                                              "enzyme_deep_amplification" = set1[1], "enzyme_deep_deletion"  = set1[2],
                                              "enzyme_shallow_amplification" = "#FB9A99", "enzyme_shallow_deletion"  = "#A6CEE3",
                                              "control" = "grey"))
        p
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_phog_en_vs_pho_sub.pdf")
        ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
        
        
        tab2p$y <- as.vector(tab2p$pho_sub)
        pos <- position_jitter(width = 0.5, seed = 1)
        p = ggplot(tab2p, aes(x=x, y=y))
        p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
        p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
        p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                force = 1, segment.size = 0.5, segment.alpha = 0.2, size=2.5,alpha=0.8, position = pos)
        p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_deep_amplification", "enzyme_deep_deletion", "enzyme_shallow_amplification", "enzyme_shallow_deletion", "control"),
                                 label = c(paste0(enzyme, "_mutation"),
                                           paste0(enzyme, "_deep_amplification"), paste0(enzyme, "_deep_deletion"),
                                           paste0(enzyme, "_shallow_amplification"), paste0(enzyme, "_shallow_deletion"),
                                           "control"))
        p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio) in", cancer))
        p = p + theme_nogrid()
        p = p + theme(axis.title.x = element_blank(),
                      axis.text.x = element_text(size= 10, vjust=0.5, hjust = 0.5, angle = 90),
                      axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
        p = p + theme(title = element_text(size = 8))
        p = p + scale_color_manual(values = c("enzyme_mutation" = "purple",
                                              "enzyme_deep_amplification" = set1[1], "enzyme_deep_deletion"  = set1[2],
                                              "enzyme_shallow_amplification" = "#FB9A99", "enzyme_shallow_deletion"  = "#A6CEE3",
                                              "control" = "grey"))
        p = p + stat_compare_means(comparisons = my_comparisons[sample_type2test], label.y = c(quantile(x = tab2p$y, probs = 0.99, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.9, na.rm = T), quantile(x = tab2p$y, probs = 0.1, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.75, na.rm = T), quantile(x = tab2p$y, probs = 0.25, na.rm = T)))
        p = p + stat_compare_means(label.y = max(tab2p$y, na.rm = T))     # Add global Anova p-value
        p
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pho_sub.pdf")
        ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
        
        
        tab2p$y <- as.vector(tab2p$pro_en)
        pos <- position_jitter(width = 0.5, seed = 1)
        p = ggplot(tab2p, aes(x=x, y=y))
        p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
        p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
        p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                force = 1, segment.size = 0.5, segment.alpha = 0.2, size=2.5,alpha=0.8, position = pos)
        p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_deep_amplification", "enzyme_deep_deletion", "enzyme_shallow_amplification", "enzyme_shallow_deletion", "control"),
                                 label = c(paste0(enzyme, "_mutation"),
                                           paste0(enzyme, "_deep_amplification"), paste0(enzyme, "_deep_deletion"),
                                           paste0(enzyme, "_shallow_amplification"), paste0(enzyme, "_shallow_deletion"),
                                           "control"))
        p = p + labs(y=paste0(enzyme, " protein abundance(log2 ratio) in", cancer))
        p = p + theme_nogrid()
        p = p + theme(axis.title.x = element_blank(),
                      axis.text.x = element_text(size= 10, vjust=0.5, hjust = 0.5, angle = 90),
                      axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
        p = p + theme(title = element_text(size = 8))
        p = p + scale_color_manual(values = c("enzyme_mutation" = "purple",
                                              "enzyme_deep_amplification" = set1[1], "enzyme_deep_deletion"  = set1[2],
                                              "enzyme_shallow_amplification" = "#FB9A99", "enzyme_shallow_deletion"  = "#A6CEE3",
                                              "control" = "grey"))
        p = p + stat_compare_means(comparisons = my_comparisons[sample_type2test], label.y = c(quantile(x = tab2p$y, probs = 0.99, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.9, na.rm = T), quantile(x = tab2p$y, probs = 0.1, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.75, na.rm = T), quantile(x = tab2p$y, probs = 0.25, na.rm = T)))
        p = p + stat_compare_means(label.y = max(tab2p$y, na.rm = T))     # Add global Anova p-value
        p
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_enzyme.pdf")
        ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
        
        tab2p$y <- as.vector(tab2p$pro_sub)
        pos <- position_jitter(width = 0.5, seed = 1)
        p = ggplot(tab2p, aes(x=x, y=y))
        p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
        p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
        p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                force = 1, segment.size = 0.5, segment.alpha = 0.2, size=2.5,alpha=0.8, position = pos)
        p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_deep_amplification", "enzyme_deep_deletion", "enzyme_shallow_amplification", "enzyme_shallow_deletion", "control"),
                                 label = c(paste0(enzyme, "_mutation"),
                                           paste0(enzyme, "_deep_amplification"), paste0(enzyme, "_deep_deletion"),
                                           paste0(enzyme, "_shallow_amplification"), paste0(enzyme, "_shallow_deletion"),
                                           "control"))
        p = p + labs(y=paste0(substrate, " protein abundance(log2 ratio) in", cancer))
        p = p + theme_nogrid()
        p = p + theme(axis.title.x = element_blank(),
                      axis.text.x = element_text(size= 10, vjust=0.5, hjust = 0.5, angle = 90),
                      axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
        p = p + theme(title = element_text(size = 8))
        p = p + scale_color_manual(values = c("enzyme_mutation" = "purple",
                                              "enzyme_deep_amplification" = set1[1], "enzyme_deep_deletion"  = set1[2],
                                              "enzyme_shallow_amplification" = "#FB9A99", "enzyme_shallow_deletion"  = "#A6CEE3",
                                              "control" = "grey"))
        p = p + stat_compare_means(comparisons = my_comparisons[sample_type2test], label.y = c(quantile(x = tab2p$y, probs = 0.99, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.9, na.rm = T), quantile(x = tab2p$y, probs = 0.1, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.75, na.rm = T), quantile(x = tab2p$y, probs = 0.25, na.rm = T)))
        p = p + stat_compare_means(label.y = max(tab2p$y, na.rm = T))     # Add global Anova p-value
        p
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_sub.pdf")
        ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
        
        tab2p$y <- as.vector(tab2p$phog_en)
        pos <- position_jitter(width = 0.5, seed = 1)
        p = ggplot(tab2p, aes(x=x, y=y))
        p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
        p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
        p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                force = 1, segment.size = 0.5, segment.alpha = 0.2, size=2.5,alpha=0.8, position = pos)
        p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_deep_amplification", "enzyme_deep_deletion", "enzyme_shallow_amplification", "enzyme_shallow_deletion", "control"),
                                 label = c(paste0(enzyme, "_mutation"),
                                           paste0(enzyme, "_deep_amplification"), paste0(enzyme, "_deep_deletion"),
                                           paste0(enzyme, "_shallow_amplification"), paste0(enzyme, "_shallow_deletion"),
                                           "control"))
        p = p + labs(y=paste0(enzyme, " phosphoprotein abundance(log2 ratio) in", cancer))
        p = p + theme_nogrid()
        p = p + theme(axis.title.x = element_blank(),
                      axis.text.x = element_text(size= 10, vjust=0.5, hjust = 0.5, angle = 90),
                      axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
        p = p + theme(title = element_text(size = 8))
        p = p + scale_color_manual(values = c("enzyme_mutation" = "purple",
                                              "enzyme_deep_amplification" = set1[1], "enzyme_deep_deletion"  = set1[2],
                                              "enzyme_shallow_amplification" = "#FB9A99", "enzyme_shallow_deletion"  = "#A6CEE3",
                                              "control" = "grey"))
        p = p + stat_compare_means(comparisons = my_comparisons[sample_type2test], label.y = c(quantile(x = tab2p$y, probs = 0.99, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.9, na.rm = T), quantile(x = tab2p$y, probs = 0.1, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.75, na.rm = T), quantile(x = tab2p$y, probs = 0.25, na.rm = T)))
        p = p + stat_compare_means(label.y = max(tab2p$y, na.rm = T))     # Add global Anova p-value
        p
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pho_enzyme.pdf")
        ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
        
        ## calculate the deviation of substrate phosphorylation of each sample comparied to the median of controls
        tab2p$pho_sub_qt <- ecdf_fun(x = tab2p$pho_sub, perc = tab2p$pho_sub)
        tab2p$pro_sub_qt <- ecdf_fun(x = tab2p$pro_sub, perc = tab2p$pro_sub)
        tab2p$phog_en_qt <- ecdf_fun(x = tab2p$phog_en, perc = tab2p$phog_en)
        tab2p$pro_en_qt <- ecdf_fun(x = tab2p$pro_en, perc = tab2p$pro_en)
        
        write.table(x = tab2p, file = paste0(subdir5, enzyme, "_", substrate, "_", rsd, ".txt"), quote = F, sep = "\t", row.names = F)
      }
      
    }
  }
}

