# Yige Wu @ WashU 2018 Apr
## check regulated pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(ggrepel)
# devtools::install_github("tidyverse/ggplot2")
library(ggplot2)

# inputs ------------------------------------------------------------------
## input druggable gene list
clinical <- fread(input = paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180307.txt"), data.table = F)

## input enzyme_substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD_extended.csv"))

mutimpact <- fread(input = paste0(ppnD, "genoalt/tables/merge_mutation_impact/mutation_impact.txt"), data.table = F)
mutimpact_pho <- mutimpact[mutimpact$substrate_type == "phosphosite",]

# set variables -----------------------------------------------------------
mut_p <- 0.05


for (cancer in c("BRCA")) {
  ## input mutation impact table
  mutimpact_pho_can <- mutimpact_pho[mutimpact_pho$Cancer == cancer & mutimpact_pho$p_value < mut_p,]
  
  for (i in 1:nrow(mutimpact_pho_can)) {
    enzyme <- as.character(mutimpact_pho_can[i, "Mutated_Gene"])
    substrate <- as.character(mutimpact_pho_can[i, "Substrate_Gene"])
    rsd <- as.character(mutimpact_pho_can[i, "SUB_MOD_RSD"])
    
    pro_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PRO_formatted_normalized_noControl.txt",sep=""), data.table = F)
    pro_data <- pro_data[pro_data$Gene == enzyme,]
    
    ## input phospho level
    pho_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
    pho_data <- pho_data[pho_data$Gene == substrate,]
    
    ## input phospho level
    phog_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt",sep=""), data.table = F)
    phog_data <- phog_data[phog_data$Gene == enzyme,]
    
    pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
    pho_data <- pho_data[pho_head$SUB_MOD_RSD == rsd,]
    
    pro.m <- melt(pro_data)
    colnames(pro.m)[ncol(pro.m)] <- "pro_kin"
    
    pho.m <- melt(pho_data)
    colnames(pho.m)[ncol(pho.m)] <- "pho_sub"
    
    phog.m <- melt(phog_data)
    colnames(phog.m)[ncol(phog.m)] <- "pho_kin"
    
    sup_tab <- merge(pro.m[, c("variable", "pro_kin")], pho.m[, c("variable", "pho_sub")], all = T)
    sup_tab <- merge(sup_tab, phog.m[, c("variable", "pho_kin")], all = T)
    sup_tab$partID <- sampID2partID(sampleID_vector = as.vector(sup_tab$variable), sample_map = clinical)
    
    ## input maf
    maf <- loadMaf(cancer = cancer, maf_files = maf_files)
    maf <- maf[(maf$Hugo_Symbol == enzyme | maf$Hugo_Symbol == substrate),]
    
    # maf <- maf[(maf$Hugo_Symbol == enzyme | maf$Hugo_Symbol == substrate | maf$Hugo_Symbol %in% unique(ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$SUB_GENE == substrate & ptms_site_pairs_sup$SUB_MOD_RSD == rsd])),]
    maf <- maf[maf$Variant_Classification != "Silent",]
    if (nrow(maf) > 0) {
      maf$partID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
      maf$aa_change <- paste0(maf$Hugo_Symbol, ":", maf$HGVSp_Short)
      maf$is.upstream <- ifelse(maf$Hugo_Symbol == enzyme, TRUE, FALSE)
      sup_tab <- merge(sup_tab, maf[, c("partID", "Variant_Classification", "aa_change", "is.upstream")], all.x = T) 
    }
    rm(maf)
    tab2p <- sup_tab
    tab2p_sort <- NULL
    for (sampID in unique(tab2p$variable)) {
      tab_tmp <- unique(tab2p[tab2p$variable == sampID, c("partID", "variable", "pro_kin", "pho_sub", "pho_kin")])
      muts <- unique(as.vector(tab2p$aa_change[tab2p$variable == sampID]))
      muts <- muts[!is.na(muts)]
      
      if (length(muts) > 0) {
        tab_tmp$aa_change <- paste0(muts, collapse = "\n")
        tab_tmp$mutated <- T
        tab_tmp$is.upstream <- any(grepl(pattern = enzyme, x = muts))
      } else {
        tab_tmp$aa_change <- NA
        tab_tmp$mutated <- F
        tab_tmp$is.upstream <- F
      }
      tab2p_sort <- rbind(tab2p_sort, tab_tmp)
    }

    pos <- position_jitter(width = 0.5, seed = 1)
    p = ggplot(tab2p_sort, aes(x=is.upstream, y=pho_sub, color = mutated, label= as.character(aa_change)))
    p = p + geom_point(aes(shape = ifelse(!is.na(is.upstream) & (is.upstream == TRUE), "b", "a")), 
                       position = pos, stroke = 0, alpha = 0.8)
    p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
    p = p + geom_text_repel(aes(segment.color = mutated),
                            force = 1,
                            segment.size = 0.5, segment.alpha = 0.2,
                            size=1.5,alpha=0.8, position = pos)
    p = p + labs(y=paste0(substrate, " ", rsd, " phosphorylation abundance(log2 ratio"))
    p = p + theme_nogrid()
    p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                  axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                  axis.text.y = element_text(colour="black", size=8))#element_text(colour="black", size=14))
    p = p + theme(title = element_text(size = 8))
    p = p + scale_color_manual(values = c("TRUE" = "firebrick1", "FALSE" = "black"))
    p
    fn = paste0(makeOutDir(resultD = resultD), enzyme, "_", substrate, "_", rsd, "_pho_sub.pdf")
    ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
  }
  
}

library(readxl)
pks <- fread(input = "./Ding_Lab/Projects_Current/IDG/IDG_shared_data/gene_lists/All_Kinase_cleaned.txt", data.table = F)
pps <- read_xlsx(path = "./pan3can_shared_data/Phospho_databases/DEPOD/DEPOD_201410_human_phosphatases.xlsx")
for (cancer in cancers_sort) {
  phog_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt",sep=""), data.table = F)
  print(nrow(phog_data[phog_data$Gene %in% pks$gene,]))
  print(nrow(phog_data[phog_data$Gene %in% pps$`Gene symbol`,]))
  
}
