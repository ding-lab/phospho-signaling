# Yige Wu @ WashU 2018 Jul
## plot how AKT1 mutated samples have substrates up-phosphorylated

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/pan3can_analysis/pan3can_aes.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# inputs ------------------------------------------------------------------
clinical <- fread(input = paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180307.txt"), data.table = F)

# process -----------------------------------------------------------------


# trans-regulated pairs overlap with driver mutations --------------
# for (cancer in c("BRCA", "OV", "CO")) {
for (cancer in c("BRCA")) {
  ## input super table
  ksea_diffexp <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  ksea_diffexp$pair <- paste0(ksea_diffexp$GENE, ":", ksea_diffexp$SUB_GENE, ":", ksea_diffexp$SUB_MOD_RSD)
  tmp <- ksea_diffexp[ksea_diffexp$FDR_pho_kin < 0.2 & !is.na(ksea_diffexp$FDR_pho_kin) & (ksea_diffexp$enzyme.driver_type != "notdriver" | ksea_diffexp$substrate.driver_type != "notdriver"),]
  tmp <- tmp[(!is.na(tmp$p_value) & tmp$p_value < 0.2) | (!is.na(tmp$p_value.pro) & tmp$p_value.pro < 0.2),]
  
  phog_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt",sep=""), data.table = F)
  pho_data <- fread(input = paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt"), data.table = F)
  pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
  maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/CPTAC2_", cancer, "_prospective.v1.3.somatic.variants.031918.maf"), data.table = F)
  
  if (nrow(tmp) > 0) {
    kegg <- read.delim(file = './PhosphoDrug/PhosphoDrug_shared_data/KEGG/hsa04151	PI3K-Akt signaling pathway.txt', header = F, col.names = "gene")
    
    # for (enzyme in unique(tmp$GENE)) {
    for (enzyme in c("AKT1")) {
      substrates <- as.character(tmp$SUB_GENE[tmp$GENE == enzyme])
      substrates <- intersect(substrates, kegg$gene)
      # for (gene in c("MAP2K4")) {
      # for (rsd in c("S268")) {
      if (enzyme %in% kegg$gene && length(substrates) > 0 && (enzyme %in% phog_data$Gene)) {
        ## input kinase phospho level
        phog_kin <- phog_data[phog_data$Gene == enzyme,]
        
        phog.m <- melt(phog_kin)
        colnames(phog.m)[ncol(phog.m)] <- "exp_value"
        phog.m$site <- "global"
        
        sup_tab <- phog.m
        sup_tab$partID <- sampID2partID(sampleID_vector = as.vector(sup_tab$variable), sample_map = clinical)
        
        for (substrate in substrates) {
          for (rsd in unique(tmp$SUB_MOD_RSD[tmp$SUB_GENE == substrate])) {
            ## input substrate phospho level
            pho_sub <- pho_data[pho_data$Gene == substrate & pho_head$SUB_MOD_RSD == rsd,]
            
            if (nrow(pho_sub) > 0) {
              pho.m <- melt(pho_sub)
              colnames(pho.m)[ncol(pho.m)] <- "exp_value"
              pho.m$site <- rsd
              pho.m$partID <- sampID2partID(sampleID_vector = as.vector(pho.m$variable), sample_map = clinical)
              
              sup_tab <- rbind(sup_tab, pho.m[,colnames(sup_tab)])
            }
          }
        }
        maf_tmp <- maf[(maf$Hugo_Symbol %in% c(enzyme, substrates)) & maf$Variant_Classification != "Silent",]
        if (nrow(maf_tmp) > 0) {
          maf_tmp$partID <- str_split_fixed(string = maf_tmp$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
          maf_tmp$aa_change <- paste0(maf_tmp$Hugo_Symbol, ":", maf_tmp$HGVSp_Short)
          maf_tmp$is.upstream <- ifelse(maf_tmp$Hugo_Symbol == enzyme, TRUE, FALSE)
          sup_tab <- merge(sup_tab, maf_tmp[, c("partID", "Variant_Classification", "aa_change", "is.upstream", "Hugo_Symbol")], all.x = T) 
        }
        
        sup_tab$is.upstream[is.na(sup_tab$is.upstream)] <- "NA"
        sup_enz <- sup_tab[sup_tab$Gene == enzyme,]
        sup_tab$variable <- factor(sup_tab$variable, levels = as.vector(sup_enz$variable)[order(sup_enz$exp_value)])
        lim = 3
        sup_tab$exp_value_capped <- sup_tab$exp_value
        sup_tab$exp_value_capped[sup_tab$exp_value_capped > lim] <- lim
        sup_tab$exp_value_capped[sup_tab$exp_value_capped < (-lim)] <- (-lim)
        
        p <- ggplot(data = sup_tab)
        p <- p + geom_tile(mapping = aes(x = variable, y = Gene, fill = exp_value_capped))
        p = p + scale_fill_gradientn(name= "Phosphorylation", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
        p = p + facet_grid(.~aa_change, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
        p = p + theme_nogrid() + theme(strip.text.x = element_text(size = 8, angle = 90), axis.text.x = element_blank())
        p
        fn = paste0(makeOutDir(resultD = resultD), enzyme, "_pro_sub+pho_kin~pho_sub.pdf")
        ggsave(file=fn, height=3, width=6, useDingbats=FALSE)
      }
    }
  }
  
}
