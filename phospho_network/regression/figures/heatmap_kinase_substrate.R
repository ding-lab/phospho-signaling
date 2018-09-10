source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(dplyr)
library(ggrepel)

# variables ---------------------------------------------------------------
reg_nonNA2test <- c(25)
top_kinase2show <- 10
top_substrate2show <- 2
# inputs -------------------------------------------------------------------
## input driver genes
driver_genes <- read.delim(file = "./TCGA_data/reference_files/Consensus.Genelist.full.txt")

## annotate kinase substrate regulation
sup_cans_tab_en$regulated <- (sup_cans_tab_en$coef_sig & sup_cans_tab_en$fdr_sig)
for (reg_nonNA in reg_nonNA2test) {
  for (enzyme_type in c("kinase")) {
    sup_cans_tab_en <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                            enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                                            "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    
    sup_cans_tab_en_reg <- sup_cans_tab_en[sup_cans_tab_en$regulated,]
    for (cancer in "OV") {
      for (enzyme in "IKBKB") {
        for (substrate in "TP53") {
          for (rsd in "S315") {
            rsd <- "S315"
            pho <- loadPhosphositeNormalizedTumor(cancer = cancer)
            pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho$Phosphosite), gene_vector = as.vector(pho$Gene))
            pho_sub <- pho[pho_head$SUBSTRATE == substrate & pho_head$SUB_MOD_RSD == rsd,]
            
            sup_tab <- melt(pho_sub)
            colnames(sup_tab) <- c("Gene", "Phosphosite", "variable", "value")
            sup_tab$x_print <- paste0(substrate, "_", rsd)
            
            pho_collapsed <- loadPhosphoproteinNormalizedTumor(cancer = cancer)
            pho_en <- pho_collapsed[pho_collapsed$Gene == enzyme,]
            pho_en.m <- melt(pho_en)
            colnames(pho_en.m) <- c("enzyme", "variable", "value")
            pho_en.m$x_print <- paste0(enzyme, "_phosphoprotein")
            sup_tab <- rbind(sup_tab[, c("variable", "value", "x_print")], pho_en.m[, c("variable", "value", "x_print")])
            sup_tab$variable <- factor(sup_tab$variable, levels = as.vector(pho_en.m$variable[order(pho_en.m$value)]))
            
            maf <- loadMaf(cancer = cancer, maf_files = maf_files)
            maf <- maf[maf$Hugo_Symbol == substrate & maf$Variant_Classification != "Silent",]
            maf$partID <- str_split_fixed(string = as.vector(maf$Tumor_Sample_Barcode), pattern = "_", 2)[,1]
            sup_tab$partID <- sampID2partID(sampleID_vector = sup_tab$variable, sample_map = loadSampMap())
            sup_tab <- merge(sup_tab, maf[, c("partID", "Hugo_Symbol","Variant_Classification")], all.x = T)
            sup_tab$x_print <- factor(sup_tab$x_print, levels = c(paste0(substrate, "_", rsd), paste0(enzyme, "_phosphoprotein")))
            p = ggplot(sup_tab)
            p = p + geom_tile(aes(x=variable, y=x_print, fill=value), width = 0.7, height = 0.7, size=0.5)#, linetype="blank")
            p = p + scale_fill_gradientn(name= "phospho_level", na.value=NA, colours=RdBu1024, limit=c(-3,3))
            p = p + theme_bw() + theme_nogrid()
            p = p + theme(axis.text.x = element_blank())
            p = p + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7))
            p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
            p
            ggsave(filename = paste0(makeOutDir(resultD = resultD), enzyme, "_", substrate, "_", rsd, ".pdf"), 
                   width = 10, height = 1.5)
          }
        }
      }
      
    }
  }
}