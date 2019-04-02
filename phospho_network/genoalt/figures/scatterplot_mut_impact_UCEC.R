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

my_comparisons <- list( enzyme_mutation = c("enzyme_mutation", "control"))
my_comparisons <- list( enzyme_mutation = c("enzyme_mutation", "control"), 
                        enzyme_mutation2 = c("enzyme_mutation", "normal"), 
                        control = c("control", "normal"))
sample_type2test <- names(my_comparisons)

# inputs ------------------------------------------------------------------

clinical <- fread(input = paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180307.txt"), data.table = F)
msi_score <- fread(input = "./cptac2p/analysis_results/phospho_network/genoalt/figures/scatterplot_mutimpact_example/CPTAC.MSI.score.tsv", data.table = F)
partIDs_msih <- msi_score$Sample[msi_score$Score >= 3.5]
partIDs_msil <- msi_score$Sample[msi_score$Score < 3.5 & msi_score$Score >= 1.0]
partIDs_mss <- msi_score$Sample[msi_score$Score < 1.0]

## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= networkin_thres),]
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


# set variables -----------------------------------------------------------
sample_type <- "tumor"


for (cancer in cancers2process) {
  ## see what's not covered in the consensus gene list
  print(SMGs[[cancer]][!(SMGs[[cancer]] %in% loadGeneList(gene_type = "driver", cancer = cancer, is.soft.limit = ""))])
  
  ## see how many SMGs per cancer types
  print(length(SMGs[[cancer]]))
  
  ## see how many SMGs are kinases
  print(length(SMGs[[cancer]][SMGs[[cancer]] %in% all_kinases$gene]))
  print(SMGs[[cancer]][SMGs[[cancer]] %in% all_kinases$gene])
}

# mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_cna_impact_UCEC/mut_cnv_sig_cans.txt"), data.table = F)
mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_complex_UCEC/mut_cnv_sig_cans.txt"), data.table = F)
mut_cnv_cans$SUB_MOD_RSD <- str_split_fixed(string = mut_cnv_cans$pair, pattern = "\\:", n = 3)[,3]

## input the complete table of mutation impact
mut_cnv_cans$pair_pro <- paste0(mut_cnv_cans$GENE, ":", mut_cnv_cans$SUB_GENE)
mut_cnv_cans$pair <- paste0(mut_cnv_cans$GENE, ":", mut_cnv_cans$SUB_GENE, ":", mut_cnv_cans$SUB_MOD_RSD)
# mut_cnv_cans <- merge(mut_cnv_cans, unique(ptms_site_pairs_sup[, c("GENE", "enzyme_type")]), all.x = T)

##filter out bad annotated enzyme type
# mut_cnv_cans <- mut_cnv_cans[!(mut_cnv_cans$GENE == "PTPRG" & mut_cnv_cans$enzyme_type == "kinase"),]

# plot mutation impact -----------------------------------------------------------
resultDnow <- makeOutDir(resultD = resultD)
subdir1 <- paste0(resultDnow, "mutation", "/")
dir.create(subdir1)

for (cancer in "UCEC") {
  subdir2 <- paste0(subdir1, cancer, "/")
  dir.create(subdir2)
  
  # maf <- loadMaf(cancer = cancer, maf_files = maf_files)
  maf <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/Baylor_DataFreeze_V2/UCEC_somatic_mutation_site_level_V2.0.maf", data.table = F)
  
  ## input somatic mutation matrix
  # mut_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_mutation_matrix_UCEC/", cancer, "_somatic_mutation_in_enzyme_substrate.txt"), data.table = F)
  mut_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_mutation_matrix_for_signor_complex_UCEC/", cancer, "_somatic_mutation_in_complex.txt"), data.table = F)
  
  mut_cnv_cans <- load_mut_impact_proteome()
  mut_cnv_tab <- mut_cnv_cans[mut_cnv_cans$p < 0.05 & mut_cnv_cans$cancer == cancer,]
  
  # mut_cnv_tab <- mut_cnv_cans
  # tab_mut <- mut_cnv_tab[mut_cnv_tab$SUB_GENE %in% c("BAD", "GSK3B", "PIK3R1", "CTNNB1", "MAP2K1", "MAPK1", "AKT1", "APC") & mut_cnv_tab$GENE %in% c("AKT1", "PIK3CA", "APC", "MAP3K1", "CTNNB1"),]
  
  for (self in c( "trans")) {
    # tab_mut <- mut_cnv_tab[mut_cnv_tab$geneA %in% c("CTNNB1") & mut_cnv_tab$geneB %in% c("APC"),]
    # tab_mut <- mut_cnv_tab[mut_cnv_tab$SUB_GENE %in% c("TP53BP1") & mut_cnv_tab$GENE %in% c("TP53") & mut_cnv_tab$SELF == self,]
    # tab_mut <- mut_cnv_tab[mut_cnv_tab$geneB %in% c("TP53BP1") & mut_cnv_tab$geneA %in% c("TP53"),]
    # tab_mut <- mut_cnv_tab[(mut_cnv_tab$geneA == "TP53" & mut_cnv_tab$geneB == "TP53BP1" & mut_cnv_tab$SUB_MOD_RSD == "S1763") | (mut_cnv_tab$geneA == "CTNNB1" & mut_cnv_tab$geneB == "APC" & mut_cnv_tab$SUB_MOD_RSD == "S2106") | (mut_cnv_tab$geneA == "CTNNB1" & mut_cnv_tab$geneB == "APC" & mut_cnv_tab$SUB_MOD_RSD == "S2278"),]
    tab_mut <- mut_cnv_tab[mut_cnv_tab$GENE == "CTNNB1" & mut_cnv_tab$SUB_GENE == "AXIN1" & mut_cnv_tab$SUB_MOD_RSD %in% c("S77", "S493"),]
    # tab_mut <- mut_cnv_tab[mut_cnv_tab$GENE == "CTNNB1" & mut_cnv_tab$SUB_GENE == "APC",]
    
    ## input protein and phosphorylation data
    # pro_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_proteomics_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
    # pho_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_phosphoproteomics_site_level_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
    # phog_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_phosphoproteomics_gene_level_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
    
    pro_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/Baylor_DataFreeze_V2/UCEC_proteomics_PNNL_ratio_median_polishing_log2_V2.0.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
    pho_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/Baylor_DataFreeze_V2/UCEC_phosphoproteomics_PNNL_ratio_median_polishing_site_level_log2_V2.0.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
    phog_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/Baylor_DataFreeze_V2/UCEC_phosphoproteomics_PNNL_ratio_median_polishing_gene_level_log2_V2.0.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
    pro_data <- data.frame(pro_data)
    colnames(pro_data)[1] <- "Gene"
    pho_data <- data.frame(pho_data)
    phog_data <- data.frame(phog_data)
    colnames(phog_data)[1] <- "Gene"

    pho_head <- data.frame(str_split_fixed(string = pho_data$idx, pattern = "-", n = 2))
    colnames(pho_head) <- c("SUBSTRATE", "SUB_MOD_RSD")
    
    ## input meta table and transform the ids
    meta_tab <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/Baylor_DataFreeze_V2/UCEC_CPTAC3_meta_table_v2.0.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)
    meta_tab <- data.frame(meta_tab)
    rownames(meta_tab) <- as.vector(meta_tab$idx)
    
    samples <- colnames(pro_data)[!(colnames(pro_data) %in% c("Gene")) & !(colnames(pro_data) %in% c("idx")) & !grepl(pattern = "Mix", x = (colnames(pro_data)))  & !grepl(pattern = "rep", x = (colnames(pro_data)))]
    samples <- samples[!(grepl(pattern = "rep", x = meta_tab[samples, "Proteomics_Participant_ID"])) & !(meta_tab[samples, "Proteomics_Participant_ID"] %in% c("C3L-00084", "C3L-01284", "C3N-01001"))]
    # samples_T <- samples[grepl(pattern = paste0('\\.', toupper(substr(x = "tumor", start = 1, stop = 1))), x = samples)]
    # samples_N <- samples[grepl(pattern = paste0('\\.', toupper(substr(x = "normal", start = 1, stop = 1))), x = samples)]
    samples_T <- samples[meta_tab[samples, "Proteomics_Tumor_Normal"] == "Tumor"]
    samples_N <- samples[meta_tab[samples, "Proteomics_Tumor_Normal"] == "Adjacent_normal"]
    samples <- c(samples_T, samples_N)
    sampIDs <- samples
    # tmp <- str_split_fixed(string = samples, pattern = "\\.", n = 3)
    
    # partIDs <- paste0(tmp[,1], "-", tmp[,2])
    partIDs <- as.vector(meta_tab[sampIDs, "Proteomics_Participant_ID"])
    names(partIDs) <- samples
    print(length(partIDs))
    # pho_data <- pho_data[, sampIDs]
    # colnames(pho_data) <- partIDs
    
    if (nrow(tab_mut) == 0) {
      print(paste0(cancer, "no data!"))
    } else {
      for (i in 1:nrow(tab_mut)) {
        enzyme <- as.character(tab_mut[i, "GENE"])
        substrate <- as.character(tab_mut[i, "SUB_GENE"])
        rsd <- as.character(tab_mut[i, "SUB_MOD_RSD"])
        
        # enzyme <- as.character(tab_mut[i, "geneA"])
        # substrate <- as.character(tab_mut[i, "geneB"])
        # rsd <- as.character(tab_mut[i, "SUB_MOD_RSD"])
        
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
        sup_tab$partID <- partIDs[as.vector(sup_tab$variable)]
        
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
          maf_en$aa_change <- paste0(maf_en$Hugo_Symbol, "\n", maf_en$HGVSp_Short)
          maf_en$is.upstream <- ifelse(maf_en$Hugo_Symbol == enzyme, TRUE, FALSE)
          maf_en$position <- str_split_fixed(string = str_split_fixed(string = maf_en$HGVSp_Short, pattern = 'p.[A-Z]', 2)[,2], pattern = '\\*|[A-Z]', n = 2)[,1]
          maf_en$position <- as.numeric(as.vector(maf_en$position))
          sup_tab <- merge(sup_tab, maf_en[, c("partID", "Variant_Classification", "aa_change", "is.upstream", "position")], all.x = T)
        }
        if (self == "trans") {
          maf_sub <- maf[(maf$Hugo_Symbol == substrate),];  maf_sub <- maf_sub[maf_sub$Variant_Classification != "Silent",]
          if (nrow(maf_sub) > 0) {
            maf_sub$partID <- str_split_fixed(string = maf_sub$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
            maf_sub$aa_change <- paste0(maf_sub$Hugo_Symbol, "\n", maf_sub$HGVSp_Short)
            sup_tab <- merge(sup_tab, maf_sub[, c("partID", "Variant_Classification", "aa_change")], by = c("partID"), all.x = T, suffixes = c("", ".sub")) 
          }
        }
        
        ## annotate sample type depending on the genomic alterations
        partIDs_overlap <- intersect(unique(sup_tab$partID), colnames(mut_mat))
        mut_mat_en <- mut_mat[mut_mat$Hugo_Symbol == enzyme,]
        if (nrow(mut_mat_en) > 0){
          mut_partIDs <- partIDs_overlap[(mut_mat_en[,partIDs_overlap] != "") & (mut_mat_en[, partIDs_overlap] != "Silent")]
        } else {
          mut_partIDs <- NULL 
        }

        sup_tab$sample_type <- "control"
        sup_tab$sample_type[sup_tab$partID %in% mut_partIDs] <- "enzyme_mutation"
        sup_tab$sample_type[sup_tab$variable %in% samples_N] <- "normal"
        
        sup_tab$subtype <- ""
        
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
        
        pos <- position_jitter(width = 0.5, seed = 1)
        tab2p <- sup_tab
        tab2p <- unique(tab2p)
        tab2p$x <- tab2p$sample_type
        tab2p$x <- factor(tab2p$x, levels = c("enzyme_mutation", "control", "enzyme_deep_amplification", "enzyme_deep_deletion", "enzyme_shallow_amplification", "enzyme_shallow_deletion", "normal"))
        tab2p$text <- NA
        tab2p$text[tab2p$sample_type == "enzyme_mutation"] <- tab2p$aa_change[tab2p$sample_type == "enzyme_mutation"]
        if (self == "trans") {
          tab2p$text[!is.na(tab2p$aa_change.sub) & is.na(tab2p$text)] <- as.vector(tab2p$aa_change.sub[!is.na(tab2p$aa_change.sub) & is.na(tab2p$text)])
          tab2p$text[!is.na(tab2p$aa_change.sub) & !is.na(tab2p$text)] <- paste0(tab2p$text[!is.na(tab2p$aa_change.sub) & !is.na(tab2p$text)], "|",  tab2p$aa_change.sub[!is.na(tab2p$aa_change.sub) & !is.na(tab2p$text)])
          tab2p$text[!is.na(tab2p$aa_change.sub) & !is.na(tab2p$text)] <- sapply(X = tab2p$text[!is.na(tab2p$aa_change.sub) & !is.na(tab2p$text)], FUN = function(s) paste0(unique(unlist(strsplit(x = s, split = "\\|"))), collapse = "|"))
        }
        tab2p$text[tab2p$sample_type == "normal"] <- NA
        
        # tab2p <- sup_tab
        tab2p <- tab2p[order(tab2p$position, tab2p$partID, decreasing = T),]
        tab2p <- tab2p[!(duplicated(tab2p[, c("partID")])),]
        
        
        tab2p$y <- as.vector(tab2p$pho_sub)
        tab2p$x <- as.vector(tab2p$sample_type)
        tab2p$x <- factor(tab2p$x, levels = c("enzyme_mutation", "control", "enzyme_deep_amplification", "enzyme_deep_deletion", "enzyme_shallow_amplification", "enzyme_shallow_deletion", "normal"))
        
        pos <- position_jitter(width = 0.5, seed = 1)
        p = ggplot(tab2p, aes(x=x, y=y))
        p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.6, size = 4)
        p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
        if (enzyme %in% c("APC", "CTNNB1", "TP53", "JAK1")) {
          p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                  force = 1, segment.size = 0.5, segment.alpha = 0.2, size = 2.5,alpha=0.6, position = pos)
        } else if (enzyme == "MAP3K1") {
          p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                  force = 1, segment.size = 0.5, segment.alpha = 0.2, size = 4, alpha=0.6, position = pos)
        } else {
          p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                  force = 1, segment.size = 0.5, segment.alpha = 0.2, size = 6, alpha=0.6, position = pos)
        }

        p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_deep_amplification", "enzyme_deep_deletion", "enzyme_shallow_amplification", "enzyme_shallow_deletion", "control", "normal"),
                                 label = c(paste0(enzyme, "_mutated\n",  "tumors"),
                                           paste0(enzyme, "_deep_amplification"), paste0(enzyme, "_deep_deletion"),
                                           paste0(enzyme, "_shallow_amplification"), paste0(enzyme, "_shallow_deletion"),
                                           "control\ntumor", "adjacent\nnormal"))
        p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio)"))
        p = p + theme_nogrid()
        p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 15, face = "bold"),
                      axis.text.x = element_text(size= 20, vjust=0.5, hjust = 0.5, face = "bold"),
                      axis.text.y = element_text(colour="black", size=8))
        p = p + theme(title = element_text(size = 18, face = "bold"))
        p = p + ggtitle(label = paste0(enzyme, " mutational association with ", substrate, "_", rsd))
        p = p + scale_color_manual(values = c("enzyme_mutation" = set1[1],
                                              "enzyme_deep_amplification" = set1[1], "enzyme_deep_deletion"  = set1[2],
                                              "enzyme_shallow_amplification" = "#FB9A99", "enzyme_shallow_deletion"  = "#A6CEE3",
                                              "control" = "purple",
                                              "normal" = "grey50"))
        p = p + stat_compare_means(comparisons = my_comparisons[sample_type2test], label.y = c(quantile(x = tab2p$y, probs = 0.99, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.9, na.rm = T), quantile(x = tab2p$y, probs = 0.1, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.75, na.rm = T), quantile(x = tab2p$y, probs = 0.25, na.rm = T)))
        p
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pho_sub.pdf")
        ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
        
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pho_sub.png")
        ggsave(file=fn, height=7, width = 8, device = png())
        dev.off()
        
        pos <- position_jitter(width = 0.5, seed = 1)
        tab2p <- tab2p[tab2p$x != "normal",]
        p = ggplot(tab2p, aes(x=x, y=y))
        p = p + geom_boxplot(aes(fill = x), alpha = 0.8)
        p = p + geom_point(aes(shape = x), position = pos, stroke = 0, alpha = 1, size = 2)
        p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_deep_amplification", "enzyme_deep_deletion", "enzyme_shallow_amplification", "enzyme_shallow_deletion", "control", "normal"),
                                 label = c(paste0(enzyme, "_mutated\n",  "tumors"),
                                           paste0(enzyme, "_deep_amplification"), paste0(enzyme, "_deep_deletion"),
                                           paste0(enzyme, "_shallow_amplification"), paste0(enzyme, "_shallow_deletion"),
                                           "control\ntumor", "adjacent\nnormal"))
        p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio)"))
        p = p + theme_nogrid()
        p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 15, face = "bold"),
                      axis.text.x = element_text(size= 20, vjust=0.5, hjust = 0.5, face = "bold"),
                      axis.text.y = element_text(colour="black", size=8))
        p = p + theme(title = element_text(size = 18, face = "bold"))
        p = p + ggtitle(label = paste0(enzyme, " mutational association with ", substrate, "_", rsd))
        p = p + scale_fill_manual(values = c("enzyme_mutation" = set1[1],
                                              "enzyme_deep_amplification" = set1[1], "enzyme_deep_deletion"  = set1[2],
                                              "enzyme_shallow_amplification" = "#FB9A99", "enzyme_shallow_deletion"  = "#A6CEE3",
                                              "control" = "purple",
                                              "normal" = "grey50"))
        p = p + stat_compare_means(comparisons = my_comparisons[sample_type2test], label.y = c(quantile(x = tab2p$y, probs = 0.99, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.9, na.rm = T), quantile(x = tab2p$y, probs = 0.1, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.75, na.rm = T), quantile(x = tab2p$y, probs = 0.25, na.rm = T)))
        p
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pho_sub_box.pdf")
        ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
        
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pho_sub_box.png")
        ggsave(file=fn, height=7, width = 8, device = png())
        dev.off()
        
        tab2p$y <- as.vector(tab2p$pro_en)
        pos <- position_jitter(width = 0.5, seed = 1)
        p = ggplot(tab2p, aes(x=x, y=y))
        p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.6, size = 4)
        p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
        if (enzyme %in% c("APC", "CTNNB1", "TP53", "JAK1")) {
          p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                  force = 1, segment.size = 0.5, segment.alpha = 0.2, size = 2.5,alpha=0.6, position = pos)
        } else if (enzyme == "MAP3K1") {
          p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                  force = 1, segment.size = 0.5, segment.alpha = 0.2, size = 4, alpha=0.6, position = pos)
        } else {
          p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                  force = 1, segment.size = 0.5, segment.alpha = 0.2, size = 6, alpha=0.6, position = pos)
        }
        p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_deep_amplification", "enzyme_deep_deletion", "enzyme_shallow_amplification", "enzyme_shallow_deletion", "control", "normal"),
                                 label = c(paste0(enzyme, "_mutated\n",  "tumors"),
                                           paste0(enzyme, "_deep_amplification"), paste0(enzyme, "_deep_deletion"),
                                           paste0(enzyme, "_shallow_amplification"), paste0(enzyme, "_shallow_deletion"),
                                           "control\ntumor", "adjacent\nnormal"))
        
        p = p + labs(y=paste0(enzyme, " protein abundance(log2 ratio) in", cancer))
        p = p + theme_nogrid()
        p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 15, face = "bold"),
                      axis.text.x = element_text(size= 20, vjust=0.5, hjust = 0.5, face = "bold"),
                      axis.text.y = element_text(colour="black", size=8))
        p = p + theme(title = element_text(size = 18, face = "bold"))
        p = p + ggtitle(label = paste0(enzyme, " mutational association with ", enzyme, " ", "protein"))
        p = p + scale_color_manual(values = c("enzyme_mutation" = set1[1],
                                              "enzyme_deep_amplification" = set1[1], "enzyme_deep_deletion"  = set1[2],
                                              "enzyme_shallow_amplification" = "#FB9A99", "enzyme_shallow_deletion"  = "#A6CEE3",
                                              "control" = "purple",
                                              "normal" = "grey50"))
        p = p + stat_compare_means(comparisons = my_comparisons[sample_type2test], label.y = c(quantile(x = tab2p$y, probs = 0.99, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.9, na.rm = T), quantile(x = tab2p$y, probs = 0.1, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.75, na.rm = T), quantile(x = tab2p$y, probs = 0.25, na.rm = T)))
        p = p + stat_compare_means(label.y = max(tab2p$y, na.rm = T))     # Add global Anova p-value
        p
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_enzyme.png")
        ggsave(file=fn, height=7, width = 8, device = png())
        dev.off()
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_enzyme.pdf")
        ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
        dev.off()
        
        tab2p$y <- as.vector(tab2p$pro_sub)
        pos <- position_jitter(width = 0.5, seed = 1)
        p = ggplot(tab2p, aes(x=x, y=y))
        p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.6, size = 4)
        p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
        if (enzyme %in% c("APC", "CTNNB1", "TP53", "JAK1")) {
          p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                  force = 1, segment.size = 0.5, segment.alpha = 0.2, size = 2.5,alpha=0.6, position = pos)
        } else if (enzyme == "MAP3K1") {
          p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                  force = 1, segment.size = 0.5, segment.alpha = 0.2, size = 4, alpha=0.6, position = pos)
        } else {
          p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                  force = 1, segment.size = 0.5, segment.alpha = 0.2, size = 6, alpha=0.6, position = pos)
        }
        p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_deep_amplification", "enzyme_deep_deletion", "enzyme_shallow_amplification", "enzyme_shallow_deletion", "control", "normal"),
                                 label = c(paste0(enzyme, "_mutated\n",  "tumors"),
                                           paste0(enzyme, "_deep_amplification"), paste0(enzyme, "_deep_deletion"),
                                           paste0(enzyme, "_shallow_amplification"), paste0(enzyme, "_shallow_deletion"),
                                           "control\ntumor", "adjacent\nnormal"))
        p = p + labs(y=paste0(substrate, " protein abundance(log2 ratio) in", cancer))
        p = p + theme_nogrid()
        p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 15, face = "bold"),
                      axis.text.x = element_text(size= 20, vjust=0.5, hjust = 0.5, face = "bold"),
                      axis.text.y = element_text(colour="black", size=8))
        p = p + theme(title = element_text(size = 18, face = "bold"))
        p = p + ggtitle(label = paste0(enzyme, " mutational association with ", substrate, " ", "protein"))
        p = p + scale_color_manual(values = c("enzyme_mutation" = set1[1],
                                              "enzyme_deep_amplification" = set1[1], "enzyme_deep_deletion"  = set1[2],
                                              "enzyme_shallow_amplification" = "#FB9A99", "enzyme_shallow_deletion"  = "#A6CEE3",
                                              "control" = "purple",
                                              "normal" = "grey50"))
        p = p + stat_compare_means(comparisons = my_comparisons[sample_type2test], label.y = c(quantile(x = tab2p$y, probs = 0.99, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.9, na.rm = T), quantile(x = tab2p$y, probs = 0.1, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.75, na.rm = T), quantile(x = tab2p$y, probs = 0.25, na.rm = T)))
        p = p + stat_compare_means(label.y = max(tab2p$y, na.rm = T))     # Add global Anova p-value
        p
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_sub.png")
        ggsave(file=fn, height=7, width = 8, device = png())
        dev.off()
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_sub.pdf")
        ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
        dev.off()
        
        pos <- position_jitter(width = 0.5, seed = 1)
        tab2p <- tab2p[tab2p$x != "normal",]
        p = ggplot(tab2p, aes(x=x, y=y))
        p = p + geom_boxplot(aes(fill = x), alpha = 0.8)
        p = p + geom_point(aes(shape = x), position = pos, stroke = 0, alpha = 1, size = 2)
        p = p + labs(y=paste0(substrate, " protein abundance(log2 ratio) in", cancer))
        p = p + theme_nogrid()
        p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 15, face = "bold"),
                      axis.text.x = element_text(size= 20, vjust=0.5, hjust = 0.5, face = "bold"),
                      axis.text.y = element_text(colour="black", size=8))
        p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_deep_amplification", "enzyme_deep_deletion", "enzyme_shallow_amplification", "enzyme_shallow_deletion", "control", "normal"),
                                 label = c(paste0(enzyme, "_mutated\n",  "tumors"),
                                           paste0(enzyme, "_deep_amplification"), paste0(enzyme, "_deep_deletion"),
                                           paste0(enzyme, "_shallow_amplification"), paste0(enzyme, "_shallow_deletion"),
                                           "control\ntumor", "adjacent\nnormal"))
        p = p + theme(title = element_text(size = 18, face = "bold"))
        p = p + ggtitle(label = paste0(enzyme, " mutational association with ", substrate, " ", "protein"))
        p = p + scale_color_manual(values = c("enzyme_mutation" = set1[1],
                                              "enzyme_deep_amplification" = set1[1], "enzyme_deep_deletion"  = set1[2],
                                              "enzyme_shallow_amplification" = "#FB9A99", "enzyme_shallow_deletion"  = "#A6CEE3",
                                              "control" = "purple",
                                              "normal" = "grey50"))
        p = p + stat_compare_means(comparisons = my_comparisons[sample_type2test], label.y = c(quantile(x = tab2p$y, probs = 0.99, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.9, na.rm = T), quantile(x = tab2p$y, probs = 0.1, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.75, na.rm = T), quantile(x = tab2p$y, probs = 0.25, na.rm = T)))
        p = p + stat_compare_means(label.y = max(tab2p$y, na.rm = T))     # Add global Anova p-value
        p

        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_sub_box.pdf")
        ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
        
        
        tab2p$y <- as.vector(tab2p$phog_en)
        pos <- position_jitter(width = 0.5, seed = 1)
        p = ggplot(tab2p, aes(x=x, y=y))
        p = p + geom_point(aes(color = sample_type, shape = subtype), position = pos, stroke = 0, alpha = 0.6, size = 4)
        p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
        if (enzyme %in% c("APC", "CTNNB1", "TP53", "JAK1")) {
          p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                  force = 1, segment.size = 0.5, segment.alpha = 0.2, size = 2.5,alpha=0.6, position = pos)
        } else if (enzyme == "MAP3K1") {
          p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                  force = 1, segment.size = 0.5, segment.alpha = 0.2, size = 4, alpha=0.6, position = pos)
        } else {
          p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = sample_type),
                                  force = 1, segment.size = 0.5, segment.alpha = 0.2, size = 6, alpha=0.6, position = pos)
        }
        p = p + scale_x_discrete(breaks = c("enzyme_mutation", "enzyme_deep_amplification", "enzyme_deep_deletion", "enzyme_shallow_amplification", "enzyme_shallow_deletion", "control", "normal"),
                                 label = c(paste0(enzyme, "_mutated\n",  "tumors"),
                                           paste0(enzyme, "_deep_amplification"), paste0(enzyme, "_deep_deletion"),
                                           paste0(enzyme, "_shallow_amplification"), paste0(enzyme, "_shallow_deletion"),
                                           "control\ntumor", "adjacent\nnormal"))
        p = p + labs(y=paste0(enzyme, " phosphoprotein abundance(log2 ratio) in", cancer))
        p = p + theme_nogrid()
        p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 15, face = "bold"),
                      axis.text.x = element_text(size= 20, vjust=0.5, hjust = 0.5, face = "bold"),
                      axis.text.y = element_text(colour="black", size=8))
        p = p + theme(title = element_text(size = 18, face = "bold"))
        p = p + ggtitle(label = paste0(enzyme, " mutational association with ", substrate, " ", "phosphorylation"))
        p = p + scale_color_manual(values = c("enzyme_mutation" = set1[1],
                                              "enzyme_deep_amplification" = set1[1], "enzyme_deep_deletion"  = set1[2],
                                              "enzyme_shallow_amplification" = "#FB9A99", "enzyme_shallow_deletion"  = "#A6CEE3",
                                              "control" = "purple",
                                              "normal" = "grey50"))
        p = p + stat_compare_means(comparisons = my_comparisons[sample_type2test], label.y = c(quantile(x = tab2p$y, probs = 0.99, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.9, na.rm = T), quantile(x = tab2p$y, probs = 0.1, na.rm = T),
                                                                                               quantile(x = tab2p$y, probs = 0.75, na.rm = T), quantile(x = tab2p$y, probs = 0.25, na.rm = T)))
        p = p + stat_compare_means(label.y = max(tab2p$y, na.rm = T))     # Add global Anova p-value
        p
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pho_enzyme.png")
        ggsave(file=fn, height=7, width = 8, device = png())
        dev.off()
        
        fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pho_enzyme.pdf")
        ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
        dev.off()

        # 
        # ## calculate the deviation of substrate phosphorylation of each sample comparied to the median of controls
        # tab2p$pho_sub_qt <- ecdf_fun(x = tab2p$pho_sub, perc = tab2p$pho_sub)
        # tab2p$pro_sub_qt <- ecdf_fun(x = tab2p$pro_sub, perc = tab2p$pro_sub)
        # tab2p$phog_en_qt <- ecdf_fun(x = tab2p$phog_en, perc = tab2p$phog_en)
        # tab2p$pro_en_qt <- ecdf_fun(x = tab2p$pro_en, perc = tab2p$pro_en)
        # 
        # write.table(x = tab2p, file = paste0(subdir5, enzyme, "_", substrate, "_", rsd, ".txt"), quote = F, sep = "\t", row.names = F)
        # 
        # stop("")
      }
      
    }
  }
}



# p = ggplot(tab2p, aes(x=position, y=pro_sub))
# p = p + geom_point(aes(color = Variant_Classification, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
# p = p + theme_nogrid()
# # p <- p + guides(color = F)
# p
# fn = paste0(subdir5, "pro_sub~position.pdf")
# ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
# 
# p = ggplot(tab2p, aes(x=position, y=pro_sub))
# p = p + geom_point(aes(color = Variant_Classification, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
# p = p + theme_nogrid()
# p <- p + xlim(c(1000, 2000))
# p
# fn = paste0(subdir5, "pro_sub~position2.pdf")
# ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
# 
# p = ggplot(tab2p, aes(x=position, y=pro_sub))
# p = p + geom_point(aes(color = Variant_Classification, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
# p = p + theme_nogrid()
# p <- p + xlim(c(1250, 1600))
# p
# fn = paste0(subdir5, "pro_sub~position3.pdf")
# ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
# 
# p = ggplot(tab2p, aes(x=position, y=phog_en))
# p = p + geom_point(aes(color = Variant_Classification, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
# p = p + theme_nogrid()
# # p <- p + guides(color = F)
# p
# fn = paste0(subdir5, "phog_en~position.pdf")
# ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
# 
# p = ggplot(tab2p, aes(x=position, y=phog_en))
# p = p + geom_point(aes(color = Variant_Classification, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
# p = p + theme_nogrid()
# p <- p + xlim(c(1000, 2000))
# p
# fn = paste0(subdir5, "phog_en~position2.pdf")
# ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
# 
# p = ggplot(tab2p, aes(x=position, y=phog_en))
# p = p + geom_point(aes(color = Variant_Classification, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
# p = p + theme_nogrid()
# p <- p + xlim(c(1250, 1600))
# p
# fn = paste0(subdir5, "phog_en~position3.pdf")
# ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
# 
# p = ggplot(tab2p, aes(x=position, y=pho_sub))
# p = p + geom_point(aes(color = Variant_Classification, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
# p = p + theme_nogrid()
# # p <- p + guides(color = F)
# p
# fn = paste0(subdir5, "pho_sub~position.pdf")
# ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
# 
# p = ggplot(tab2p, aes(x=position, y=pho_sub))
# p = p + geom_point(aes(color = Variant_Classification, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
# p = p + theme_nogrid()
# p <- p + xlim(c(1000, 2000))
# p
# fn = paste0(subdir5, "pho_sub~position2.pdf")
# ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
# 
# p = ggplot(tab2p, aes(x=position, y=pho_sub))
# p = p + geom_point(aes(color = Variant_Classification, shape = subtype), position = pos, stroke = 0, alpha = 0.8)
# p = p + theme_nogrid()
# p <- p + xlim(c(1250, 1600))
# p
# fn = paste0(subdir5, "pho_sub~position3.pdf")
# ggsave(file=fn, height=7, width = 8, useDingbats=FALSE)
