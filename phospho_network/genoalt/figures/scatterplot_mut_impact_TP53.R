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

getTP53domain <- function(position) {
  if (is.na(position)) {
    return("other")
  } else if (position <= 43) {
    return("TAD1")
  } else if (position <= 63) {
    return("TAD2")
  } else if (position <= 92) {
    return("proline rich\ndomain")
  } else  if (position >= 102 & position <= 292) {
    return("DBD")
  } else if (position >= 316 & position <= 325) {
    return("NLS")
  } else if (position >= 325 & position <= 355) {
    return("tetramerization\ndomain")
  } else if (position >= 356 & position <= 393) {
    return("regulatory\ndomain")
  } else {
    return("other")
  }
}

getTP53domains <- function(positions) {
  domains <- sapply(positions, getTP53domain)
  return(domains)
}

getPNGsize <- function(unique_number_x) {
  widths <- c(8, 8, 9, 10)
  names(widths) <- as.character(4:7)
  
  heights <- c(7, 7, 7,  7)
  names(heights) <- as.character(4:7)
  
  size <- c(widths[as.character(unique_number_x)], heights[as.character(unique_number_x)])
  return(size)
}

plot_scatterplot <- function(tab2p, pos) {
  p = ggplot(tab2p, aes(x=x, y=y))
  p = p + geom_point(aes(fill = x, shape = subtype), position = pos, stroke = 0.2, alpha = 0.6, size = 4, color = "black")
  p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
  p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text)),
                          force = 2, segment.size = 0.5, segment.alpha = 0.2, size = 2.5,alpha=0.6, position = pos, color = "black")
  p = p + scale_shape_manual(values = subtype_shape_values)
  p = p + theme_nogrid()
  p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 15, face = "bold"),
                axis.text.x = element_text(size= 20, vjust=0.5, hjust = 0.5, face = "bold"),
                axis.text.y = element_text(colour="black", size=8))
  p = p + theme(title = element_text(size = 18, face = "bold"))
  return(p)
}

plot_linearscatterplot <- function(tab2p) {
  p = ggplot(tab2p, aes(x=x, y=y))
  p = p + geom_point(aes(fill = domain, shape = subtype), stroke = 0.2, alpha = 0.6, size = 4, color = "black")
  p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text), color = variant_class),
                          force = 2, segment.size = 0.5, segment.alpha = 0.2, size = 2.5,alpha=0.6)
  p = p + scale_shape_manual(values = subtype_shape_values)
  p = p + scale_color_manual(values = c("missense" = set1[1], "truncation" = set1[2], "other" = "grey50"))
  if ("cancer" %in% colnames(tab2p)) {
    p = p + facet_grid(cancer~., scales = "free_y", space = "fixed")
  }
  p = p + theme_nogrid()
  p = p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 15, face = "bold"),
                axis.text.x = element_text(size= 15, vjust=0.5, hjust = 0.5, face = "bold"),
                axis.text.y = element_text(colour="black", size=8))
  p = p + theme(title = element_text(size = 18, face = "bold"))
  return(p)
}

# inputs ------------------------------------------------------------------
cancer_full_names <- c("breast", "ovarian", "colorectal"); names(cancer_full_names) <- cancers_sort
clinical <- fread(input = paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180307.txt"), data.table = F)
msi_score <- fread(input = "./cptac2p/analysis_results/phospho_network/genoalt/figures/scatterplot_mutimpact_example/CPTAC.MSI.score.tsv", data.table = F)
partIDs_msih <- msi_score$Sample[msi_score$Score >= 3.5]
partIDs_msil <- msi_score$Sample[msi_score$Score < 3.5 & msi_score$Score >= 1.0]
partIDs_mss <- msi_score$Sample[msi_score$Score < 1.0]

# input and parse mutaional impact result table ---------------------------
mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/parse_mut_cnv_impact/mut_cnv_cans.txt"), data.table = F)

## input the complete table of mutation impact
mut_cnv_cans$pair_pro <- paste0(mut_cnv_cans$GENE, ":", mut_cnv_cans$SUB_GENE)
mut_cnv_cans$pair <- paste0(mut_cnv_cans$GENE, ":", mut_cnv_cans$SUB_GENE, ":", mut_cnv_cans$SUB_MOD_RSD)

## input known kinase/phosphatase relationships
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)
ptms_site_pairs_sup$pair_pro_rev <- paste0(ptms_site_pairs_sup$SUB_GENE, ":", ptms_site_pairs_sup$GENE)

# plot mutation impact -----------------------------------------------------------
resultDnow <- makeOutDir(resultD = resultD)
subdir1 <- paste0(resultDnow, "mutation", "/")
dir.create(subdir1)


test <- mut_cnv_cans[mut_cnv_cans$p < 0.05 & mut_cnv_cans$GENE == "TP53",]
test <- test[order(test$p),]
# pairs2plot <- unique(mut_cnv_cans$pair[mut_cnv_cans$p < 0.05 & mut_cnv_cans$GENE == "TP53"])
# pairs2plot <- c("TP53:CDK1:Y15")
# pairs2plot <- c("TP53:PPM1D:S40", "TP53:PPM1D:T34", "TP53:PPM1D:T529")]
# pairs2plot <- unique(test$pair[test$SUB_GENE %in% c("BAX", "CDKN1A", "TIGAR", "SIRT1", "KAT5", "CREBBP", "ATM", "ATR", "CHEK1", "CHEK2")])
# pairs2plot <- unique(test$pair[test$SUB_GENE %in% c("BAX", "CDKN1A", "TIGAR")])
# pairs2plot <- unique(test$pair[test$SUB_GENE %in% c("ATM", "ATR", "CHEK1", "CHEK2")])
pairs2plot <- unique(mut_cnv_cans$pair[mut_cnv_cans$GENE == "TP53" &  mut_cnv_cans$p < 0.3 & mut_cnv_cans$SUB_GENE %in% c("PPM1D", "TP53BP1", "MDM2", "MDM4")])
pairs2plot <- unique(mut_cnv_cans$pair[mut_cnv_cans$GENE == "TP53" & mut_cnv_cans$SUB_GENE %in% c("MDM2", "MDM4")])
pairs2plot <- unique(mut_cnv_cans$pair[mut_cnv_cans$GENE == "TP53" & mut_cnv_cans$p < 0.05 &  mut_cnv_cans$SUB_GENE %in% c("TP53BP1")])

for (pair in pairs2plot) {
  enzyme <- str_split(string = pair, pattern = "\\:")[[1]][1]
  substrate <- str_split(string = pair, pattern = "\\:")[[1]][2]
  rsd <- str_split(string = pair, pattern = "\\:")[[1]][3]
  
  self <- ifelse(enzyme == substrate, "cis", "trans")
  sup_tab_cans <- NULL
  for (cancer in c("UCEC", cancers_sort)) {
    subdir2 <- paste0(subdir1, cancer, "/")
    dir.create(subdir2)
    subdir3 <- paste0(subdir2, enzyme, "/")
    dir.create(subdir3)
    subdir4 <- paste0(subdir3, substrate, "/")
    dir.create(subdir4)
    subdir5 <- paste0(subdir4, rsd, "/")
    dir.create(subdir5)
    
    ## input somatic mutation matrix
    if (cancer %in% cancers_sort) {
      mut_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_mutation_matrix/", cancer, "_somatic_mutation_in_enzyme_substrate.txt"), data.table = F)
      pro_data <- loadProteinNormalizedTumor(cancer = cancer)
      pho_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
      phog_data <- loadPhosphoproteinNormalizedTumor(cancer = cancer)
      pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
      
    } else if (cancer == "UCEC") {
      mut_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_mutation_matrix_UCEC/", cancer, "_somatic_mutation_in_enzyme_substrate.txt"), data.table = F)
      ## input protein and phosphorylation data
      pro_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_proteomics_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
      pro_data <- data.frame(pro_data)
      colnames(pro_data)[1] <- "Gene"
      
      pho_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_phosphoproteomics_site_level_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
      pho_data <- data.frame(pho_data)
      
      phog_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_phosphoproteomics_gene_level_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
      phog_data <- data.frame(phog_data)
      colnames(phog_data)[1] <- "Gene"
      
      pho_head <- data.frame(str_split_fixed(string = pho_data$idx, pattern = "-", n = 2))
      colnames(pho_head) <- c("SUBSTRATE", "SUB_MOD_RSD")
      
      samples <- colnames(pro_data)[!(colnames(pro_data) %in% c("idx")) & !grepl(pattern = "Mix", x = (colnames(pro_data)))  & !grepl(pattern = "rep", x = (colnames(pro_data)))]
      samples_T <- samples[grepl(pattern = paste0('\\.', toupper(substr(x = "tumor", start = 1, stop = 1))), x = samples)]
      samples_N <- samples[grepl(pattern = paste0('\\.', toupper(substr(x = "normal", start = 1, stop = 1))), x = samples)]
      samples <- samples_T
      sampIDs <- samples_T
      tmp <- str_split_fixed(string = samples, pattern = "\\.", n = 3)
      partIDs <- paste0(tmp[,1], "-", tmp[,2])
      names(partIDs) <- samples
      
      pho_data <- pho_data[, sampIDs]
      pro_data <- pro_data[, c("Gene", sampIDs)]
      phog_data <- phog_data[, c("Gene", sampIDs)]
    }
    maf <- loadMaf(cancer = cancer, maf_files = maf_files)
    
    ## input phospho level
    pho_sub <- pho_data[pho_head$SUBSTRATE == substrate & pho_head$SUB_MOD_RSD == rsd,]
    pho.m <- melt(pho_sub)
    colnames(pho.m)[ncol(pho.m)] <- "pho_sub"
    pho.m <- pho.m[!is.na(pho.m$pho_sub),]
    
    ## input protein level of the enzyme
    pro_en <- pro_data[pro_data$Gene == enzyme,]
    pro_en.m <- melt(pro_en)
    colnames(pro_en.m) <- c("enzyme", "variable", "pro_en")
    pro_en.m <- pro_en.m[!is.na(pro_en.m$pro_en),]
    
    ## input protein level of the substrate
    pro_sub <- pro_data[pro_data$Gene == substrate,]
    pro_sub.m <- melt(pro_sub)
    colnames(pro_sub.m) <- c("substrate", "variable", "pro_sub")
    
    sup_tab <- pho.m[, c("variable", "pho_sub")]
    sup_tab <- merge(sup_tab, pro_sub.m[, c("variable", "pro_sub")], all = T)
    sup_tab <- merge(sup_tab, pro_en.m[, c("variable", "pro_en")], all = T)
    
    sup_tab <- sup_tab[!is.na(sup_tab$pho_sub) | !is.na(sup_tab$pro_en) | !is.na(sup_tab$pro_sub),]
    if (!any(!is.na(sup_tab$pho_sub) | !is.na(sup_tab$pro_en) | !is.na(sup_tab$pro_sub))) {
      next()
    }
    
    if (cancer %in% cancers_sort) {
      sup_tab$partID <- sampID2partID(sampleID_vector = as.vector(sup_tab$variable), sample_map = clinical)
    } else {
      tmp <- str_split_fixed(string = as.vector(sup_tab$variable), pattern = "\\.", 3)
      sup_tab$partID <- paste0(tmp[,1], "-", tmp[,2])
    }
    
    ## input maf
    maf_en <- maf[(maf$Hugo_Symbol == enzyme),]
    # maf <- maf[(maf$Hugo_Symbol == enzyme | maf$Hugo_Symbol == substrate | maf$Hugo_Symbol %in% unique(ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$SUB_GENE == substrate & ptms_site_pairs_sup$SUB_MOD_RSD == rsd])),]
    maf_en <- maf_en[maf_en$Variant_Classification != "Silent",]
    if (nrow(maf_en) > 0) {
      maf_en$partID <- str_split_fixed(string = maf_en$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
      if (self == "cis") {
        maf_en$aa_change <- paste0(maf_en$HGVSp_Short)
      } else {
        maf_en$aa_change <- paste0(maf_en$Hugo_Symbol, "\n", maf_en$HGVSp_Short)
      }
      maf_en$is.upstream <- ifelse(maf_en$Hugo_Symbol == enzyme, TRUE, FALSE)
      maf_en$position <- str_split_fixed(string = str_split_fixed(string = maf_en$HGVSp_Short, pattern = 'p.[A-Z]', 2)[,2], pattern = '\\*|[A-Z]', n = 2)[,1]
      maf_en$position <- as.numeric(as.vector(maf_en$position))
      sup_tab <- merge(sup_tab, maf_en[, c("partID", "Variant_Classification", "aa_change", "is.upstream", "position")], all.x = T)
    }
    maf_sub <- maf[(maf$Hugo_Symbol == substrate),];  maf_sub <- maf_sub[maf_sub$Variant_Classification != "Silent",]
    maf_sub$partID <- str_split_fixed(string = maf_sub$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
    
    if (nrow(maf_sub) > 0) {
      if (self == "cis") {
        maf_sub$aa_change <- paste0(maf_sub$HGVSp_Short)
      } else {
        maf_sub$aa_change <- paste0(maf_sub$Hugo_Symbol, "\n", maf_sub$HGVSp_Short)
      }
      sup_tab <- merge(sup_tab, maf_sub[, c("partID", "Variant_Classification", "aa_change")], by = c("partID"), all.x = T, suffixes = c("", ".sub")) 
    } else {
      sup_tab$Variant_Classification.sub <- NA
      sup_tab$aa_change.sub <- NA
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
    sample_type2test <- data.frame(table(unique(sup_tab[, c("sample_type", "partID")])[, "sample_type"]))
    sample_type2test <- as.vector(sample_type2test$Var1[sample_type2test$Freq >= num_genoalt_thres & sample_type2test$Var1 != "control"])
    
    if (cancer == "BRCA") {
      sup_tab$subtype <- sampID2pam50(sampleID_vector = as.vector(sup_tab$variable), pam50_map = loadPAM50Map(), sample_map = loadSampMap())
      subtype_shape_values <- c("Basal" = 21, "Her2" = 22, "LumA" = 23, "LumB" = 24, "Normal" = 25, "no_RNA-seq" = 20)
    } else if (cancer == "CO") {
      sup_tab$subtype <- "other"
      sup_tab$subtype[sup_tab$partID %in% partIDs_msih] <- "MSI-High"
      sup_tab$subtype[sup_tab$partID %in% partIDs_msil] <- "MSI-Low"
      sup_tab$subtype[sup_tab$partID %in% partIDs_mss] <- "MSS"
      subtype_shape_values <- c("MSI-High" = 21, "MSI-Low" = 22, "MSS" = 23, "other" = 25)
    } else {
      sup_tab$subtype <- "other"
      subtype_shape_values <- c("other" = 25)
    }
    
    sup_tab$text <- NA
    sup_tab$text[sup_tab$sample_type == "enzyme_mutation"] <- sup_tab$aa_change[sup_tab$sample_type == "enzyme_mutation"]
    sup_tab$text[!is.na(sup_tab$aa_change.sub) & is.na(sup_tab$text)] <- as.vector(sup_tab$aa_change.sub[!is.na(sup_tab$aa_change.sub) & is.na(sup_tab$text)])
    sup_tab$text[!is.na(sup_tab$aa_change.sub) & !is.na(sup_tab$text)] <- paste0(sup_tab$text[!is.na(sup_tab$aa_change.sub) & !is.na(sup_tab$text)], "|",  sup_tab$aa_change.sub[!is.na(sup_tab$aa_change.sub) & !is.na(sup_tab$text)])
    sup_tab$text[!is.na(sup_tab$aa_change.sub) & !is.na(sup_tab$text)] <- sapply(X = sup_tab$text[!is.na(sup_tab$aa_change.sub) & !is.na(sup_tab$text)], FUN = function(s) paste0(unique(unlist(strsplit(x = s, split = "\\|"))), collapse = "|"))
    sup_tab <- unique(sup_tab)
    sup_tab <- sup_tab[!(duplicated(sup_tab[, c("partID")])),]
    # if (self == "trans") {
    #   sup_tab$text[!is.na(sup_tab$aa_change.sub) & is.na(sup_tab$text)] <- as.vector(sup_tab$aa_change.sub[!is.na(sup_tab$aa_change.sub) & is.na(sup_tab$text)])
    #   sup_tab$text[!is.na(sup_tab$aa_change.sub) & !is.na(sup_tab$text)] <- paste0(sup_tab$text[!is.na(sup_tab$aa_change.sub) & !is.na(sup_tab$text)], "|",  sup_tab$aa_change.sub[!is.na(sup_tab$aa_change.sub) & !is.na(sup_tab$text)])
    # }
    
    ## get the domain information of TP53
    sup_tab$domain <- getTP53domains(position = as.vector(sup_tab$position))
    sup_tab$domain[sup_tab$sample_type == "control"] <- "control"
    
    ## distinguish missense and truncation
    sup_tab$variant_class <- "other"
    sup_tab$variant_class[sup_tab$Variant_Classification == "Missense_Mutation"] <- "missense"
    sup_tab$variant_class[sup_tab$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation")] <- "truncation"
    
    
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
    sup_tab$cancer <- cancer
    sup_tab_cans <- rbind(sup_tab, sup_tab_cans)
    
    pos <- position_jitter(width = 0.5, seed = 1)
    if (any(!is.na(sup_tab$pho_sub))) {
      # plot pho_sub group by domain -------------------------------------------
      tab2p <- sup_tab
      tab2p$y <- as.vector(tab2p$pho_sub)
      tab2p$x <- tab2p$domain
      tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
      tab2p$x <- factor(tab2p$x, levels = c("TAD1", "TAD2", "proline rich\ndomain", "DBD", "NLS", "tetramerization\ndomain", "regulatory\ndomain", "other", "control"))
      p <- plot_scatterplot(tab2p = tab2p, pos = pos)
      p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio)"))
      p = p + ggtitle(label = paste0(enzyme, " mutational association with ", substrate, "_", rsd))
      p
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pho_sub_by_domain.pdf")
      pngsize <- getPNGsize(unique_number_x = length(unique(tab2p$x)))
      ggsave(file=fn, width = pngsize[1], height = pngsize[2])
      
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pho_sub_by_domain.png")
      ggsave(file=fn, width = pngsize[1], height = pngsize[2], device = png())
      dev.off()
      
      # plot pho_sub linear by domain -------------------------------------------
      tab2p <- sup_tab
      tab2p$y <- as.vector(tab2p$pho_sub)
      tab2p$x <- tab2p$position
      tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
      p <- plot_linearscatterplot(tab2p = tab2p)
      p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio)"))
      p = p + ggtitle(label = paste0(enzyme, " mutational association with ", substrate, "_", rsd))
      p
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pho_sub_by_domain_linear.pdf")
      pngsize <- c(10,6)
      ggsave(file=fn, width = pngsize[1], height = pngsize[2])
      
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pho_sub_by_domain_linear.png")
      ggsave(file=fn, width = pngsize[1], height = pngsize[2], device = png())
      dev.off()
    }

    if (any(!is.na(sup_tab$pro_en))) {
      # plot pro_enzyme group by domain -------------------------------------------
      tab2p <- sup_tab
      tab2p$y <- as.vector(tab2p$pro_en)
      tab2p$x <- tab2p$domain
      tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
      tab2p$x <- factor(tab2p$x, levels = c("TAD1", "TAD2", "proline rich\ndomain", "DBD", "NLS", "tetramerization\ndomain", "regulatory\ndomain", "other", "control"))
      p <- plot_scatterplot(tab2p = tab2p, pos = pos)
      p = p + labs(y=paste0(enzyme, " protein abundance(log2 ratio)"))
      p = p + ggtitle(label = paste0(enzyme, " mutational association with ", enzyme))
      p
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_enzyme_by_domain.pdf")
      pngsize <- getPNGsize(unique_number_x = length(unique(tab2p$x)))
      ggsave(file=fn, width = pngsize[1], height = pngsize[2])
      
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_enzyme_by_domain.png")
      ggsave(file=fn, width = pngsize[1], height = pngsize[2], device = png())
      dev.off()
      
      # plot pro_enzyme linear by domain -------------------------------------------
      tab2p <- sup_tab
      tab2p$y <- as.vector(tab2p$pro_en)
      tab2p$x <- tab2p$position
      tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
      p <- plot_linearscatterplot(tab2p = tab2p)
      p = p + labs(y=paste0(enzyme, " protein abundance(log2 ratio)"))
      p = p + ggtitle(label = paste0(enzyme, " mutational association with ", enzyme))
      p
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_enzyme_by_domain_linear.pdf")
      pngsize <- c(10,6)
      ggsave(file=fn, width = pngsize[1], height = pngsize[2])
      
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_enzyme_by_domain_linear.png")
      ggsave(file=fn, width = pngsize[1], height = pngsize[2], device = png())
      dev.off()
    }
    
    if (any(!is.na(sup_tab$pro_sub))) {
      # plot pro_sub group by domain -------------------------------------------
      tab2p <- sup_tab
      tab2p$y <- as.vector(tab2p$pro_sub)
      tab2p$x <- tab2p$domain
      tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
      tab2p$x <- factor(tab2p$x, levels = c("TAD1", "TAD2", "proline rich\ndomain", "DBD", "NLS", "tetramerization\ndomain", "regulatory\ndomain", "other", "control"))
      p <- plot_scatterplot(tab2p = tab2p, pos = pos)
      p = p + labs(y=paste0(substrate, " protein abundance(log2 ratio)"))
      p = p + ggtitle(label = paste0(enzyme, " mutational association with ", substrate))
      p
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_sub_by_domain.pdf")
      pngsize <- getPNGsize(unique_number_x = length(unique(tab2p$x)))
      ggsave(file=fn, width = pngsize[1], height = pngsize[2])
      
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_sub_by_domain.png")
      ggsave(file=fn, width = pngsize[1], height = pngsize[2], device = png())
      dev.off()
      
      # plot pro_enzyme linear by domain -------------------------------------------
      tab2p <- sup_tab
      tab2p$y <- as.vector(tab2p$pro_sub)
      tab2p$x <- tab2p$position
      tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
      p <- plot_linearscatterplot(tab2p = tab2p)
      p = p + labs(y=paste0(substrate, " protein abundance(log2 ratio)"))
      p = p + ggtitle(label = paste0(enzyme, " mutational association with ", substrate))
      p
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_sub_by_domain_linear.pdf")
      pngsize <- c(10,6)
      ggsave(file=fn, width = pngsize[1], height = pngsize[2])
      
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_sub_by_domain_linear.png")
      ggsave(file=fn, width = pngsize[1], height = pngsize[2], device = png())
      dev.off()
    }
    
    if (any(!is.na(sup_tab$pro_sub) & !is.na(sup_tab$pro_en))) {
      # plot pro_sub ~ pro_en -------------------------------------------
      tab2p <- sup_tab
      tab2p$y <- as.vector(tab2p$pro_sub)
      tab2p$x <- as.vector(tab2p$pro_en)
      tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
      p = ggplot(tab2p, aes(x=x, y=y))
      p = p + geom_point(aes(fill = domain, shape = subtype), stroke = 0.2, alpha = 0.6, size = 4, color = "black")
      p = p + geom_text_repel(data = tab2p, mapping = aes(segment.color = sample_type, label= as.character(text)),
                              force = 2, segment.size = 0.5, segment.alpha = 0.2, size = 2.5,alpha=0.6, color = "black")
      p = p + scale_shape_manual(values = subtype_shape_values)
      p = p + theme_nogrid()
      p = p + theme(axis.title.y = element_text(size = 15, face = "bold"),
                    axis.title.x = element_text(size = 15, face = "bold"),
                    axis.text.x = element_text(colour="black", size=8),
                    axis.text.y = element_text(colour="black", size=8))
      p = p + theme(title = element_text(size = 18, face = "bold"))
      p = p + labs(y=paste0(substrate, " protein abundance(log2 ratio)"), x=paste0(enzyme, " protein abundance(log2 ratio)"))
      p = p + ggtitle(label = paste0(enzyme, " mutational association with ", substrate, " protein"))
      p = p + guides(fill=guide_legend(override.aes=list(shape=21)))
      p
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_en~pro_sub.pdf")
      pngsize <- getPNGsize(unique_number_x = length(unique(tab2p$x)))
      ggsave(file=fn, width = 8, height = 6)
      
      fn = paste0(subdir5, enzyme, "_", substrate, "_", rsd, "_pro_en~pro_sub.png")
      ggsave(file=fn, width = 8, height = 6, device = png())
      dev.off()
    }
  }
  next()
  pos <- position_jitter(width = 0.5, seed = 1)
  pngsize <- c(10,12)
  # plot pho_sub group by domain all_cancers -------------------------------------------
  tab2p <- sup_tab_cans
  tab2p$y <- as.vector(tab2p$pho_sub)
  tab2p$x <- tab2p$domain
  tab2p$subtype <- "other"
  tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
  tab2p$x <- factor(tab2p$x, levels = c("TAD1", "TAD2", "proline rich\ndomain", "DBD", "NLS", "tetramerization\ndomain", "regulatory\ndomain", "other", "control"))
  p <- plot_scatterplot(tab2p = tab2p, pos = pos)
  p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio)"))
  p = p + ggtitle(label = paste0(enzyme, " mutational association with ", substrate, "_", rsd))
  p
  fn = paste0(subdir1, enzyme, "_", substrate, "_", rsd, "_pho_sub_by_domain.pdf")
  pngsize <- getPNGsize(unique_number_x = length(unique(tab2p$x)))
  ggsave(file=fn, width = pngsize[1], height = pngsize[2])
  
  fn = paste0(subdir1, enzyme, "_", substrate, "_", rsd, "_pho_sub_by_domain.png")
  ggsave(file=fn, width = pngsize[1], height = pngsize[2], device = png())
  dev.off()
  
  # plot pro_enzyme group by domain all cancers -------------------------------------------
  # tab2p <- sup_tab_cans
  # tab2p$y <- as.vector(tab2p$pro_en)
  # tab2p$x <- tab2p$domain
  # tab2p$subtype <- "other"
  # tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
  # tab2p$x <- factor(tab2p$x, levels = c("TAD1", "TAD2", "proline rich\ndomain", "DBD", "NLS", "tetramerization\ndomain", "regulatory\ndomain", "other", "control"))
  # p <- plot_scatterplot(tab2p = tab2p, pos = pos)
  # p = p + labs(y=paste0(enzyme, " protein abundance(log2 ratio)"))
  # p = p + ggtitle(label = paste0(enzyme, " mutational association with ", enzyme))
  # p
  # fn = paste0(subdir1, enzyme, "_", substrate, "_", rsd, "_pro_enzyme_by_domain.pdf")
  # pngsize <- getPNGsize(unique_number_x = length(unique(tab2p$x)))
  # ggsave(file=fn, width = pngsize[1], height = pngsize[2])
  # 
  # fn = paste0(subdir1, enzyme, "_", substrate, "_", rsd, "_pro_enzyme_by_domain.png")
  # ggsave(file=fn, width = pngsize[1], height = pngsize[2], device = png())
  # dev.off()
  
  # plot pro_enzyme linear by domain all cancers -------------------------------------------
  # tab2p <- sup_tab_cans
  # tab2p$y <- as.vector(tab2p$pro_en)
  # tab2p$x <- tab2p$position
  # tab2p$subtype <- "other"
  # tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
  # p <- plot_linearscatterplot(tab2p = tab2p)
  # p = p + labs(y=paste0(enzyme, " protein abundance(log2 ratio)"))
  # p = p + ggtitle(label = paste0(enzyme, " mutational association with ", enzyme, " protein"))
  # p
  # fn = paste0(subdir1, enzyme, "_", substrate, "_", rsd, "_pro_enzyme_by_domain_linear.pdf")
  # ggsave(file=fn, width = pngsize[1], height = pngsize[2])
  # 
  # fn = paste0(subdir1, enzyme, "_", substrate, "_", rsd, "_pro_enzyme_by_domain_linear.png")
  # ggsave(file=fn, width = pngsize[1], height = pngsize[2], device = png())
  # dev.off()
  
  # plot pho_sub linear by domain all cancers-------------------------------------------
  tab2p <- sup_tab_cans
  tab2p$y <- as.vector(tab2p$pho_sub)
  tab2p$x <- tab2p$position
  tab2p$subtype <- "other"
  tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
  p <- plot_linearscatterplot(tab2p = tab2p)
  p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio)"))
  p = p + ggtitle(label = paste0(enzyme, " mutational association with ", substrate, "_", rsd))
  p
  fn = paste0(subdir1, enzyme, "_", substrate, "_", rsd, "_pho_sub_by_domain_linear.pdf")
  ggsave(file=fn, width = pngsize[1], height = pngsize[2])
  
  fn = paste0(subdir1, enzyme, "_", substrate, "_", rsd, "_pho_sub_by_domain_linear.png")
  ggsave(file=fn, width = pngsize[1], height = pngsize[2], device = png())
  dev.off()
}
