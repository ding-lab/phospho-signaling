# Yige Wu @ WashU 2018 Aug
# barplot for the top kinases with high and low kinase-substrate pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R')

library(dplyr)
library(ggrepel)
library(biomaRt)
# variables ---------------------------------------------------------------
reg_nonNA2test <- c(25)
top_kinase2show <- 10
top_substrate2show <- 2
pdf_sizes <- list(trans = list(BRCA = list(width = 7, height = 6),
                               OV = list(width = 4, height = 3),
                               CO = list(width = 5, height = 4)))
pdf_sizes_per_yfacet <- list("1" = 2.5, "2" = 4, "3" = 5, "4" = 6, "5" = 7, "16" = 20, "20" = 25)

# inputs -------------------------------------------------------------------
## input driver genes
# driver_genes <- read.delim(file = "./TCGA_data/reference_files/Consensus.Genelist.full.txt")
pho_head_cancers_sim <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid.txt"), data.table = F)
pho_head_cancers_sim <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid_hg38_hg19.txt"), data.table = F)

driver_pancan <- loadGeneList(gene_type = "driver", cancer = "PANCAN", is.soft.limit = "")

ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org") ## 38
filters = listFilters(ensembl)

pam50_map <- loadPAM50Map()
## annotate kinase substrate regulation
sup_cans_tab_en$regulated <- (sup_cans_tab_en$coef_sig & sup_cans_tab_en$fdr_sig)
for (reg_nonNA in reg_nonNA2test) {
  for (enzyme_type in c("kinase")) {
    sup_cans_tab_en <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                            enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                                            "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    
    sup_cans_tab_en_reg <- sup_cans_tab_en[sup_cans_tab_en$regulated,]
    for (self in "trans") {
      for (cancer in "BRCA") {
        for (sub_gene in "ERBB2") {
          tab_gene <- sup_cans_tab_en_reg
          tab_gene <- tab_gene[tab_gene$SUB_GENE == gene & tab_gene$SELF == self & tab_gene$Cancer == cancer,]
          tab_gene <- merge(tab_gene, 
                            pho_head_cancers_sim[pho_head_cancers_sim$SUBSTRATE == gene & pho_head_cancers_sim$SUB_MOD_RSD %in% tab_gene$SUB_MOD_RSD, 
                                                 c("SUBSTRATE", "SUB_MOD_RSD", "seq_region_name", "start_hg19")], 
                            by.x = c("SUB_GENE", "SUB_MOD_RSD"), by.y = c("SUBSTRATE", "SUB_MOD_RSD"),
                            all.x = T)
          refseqs_mapTab = getBM(attributes = c("refseq_mrna", "refseq_peptide", "hgnc_symbol"),
                                 filters = c("refseq_peptide"),
                                 values = str_split_fixed(string = as.vector(tab_gene$transcript), pattern = "\\.", 2)[,1], 
                                 mart = ensembl, uniqueRows=T)
          
          
          tab4proteinpaint <- data.frame(gene = tab_gene$SUB_GENE, 
                                         refseq = refseqs_mapTab$refseq_mrna, 
                                         chromosome = tab_gene$seq_region_name, 
                                         start = tab_gene$start_hg19, 
                                         aachange = tab_gene$SUB_MOD_RSD, class = "silent")
          tab4proteinpaint <- tab4proteinpaint[!is.na(tab4proteinpaint$start),]
          tab4proteinpaint <- unique(tab4proteinpaint)
          write.table(x = tab4proteinpaint, file = paste0(makeOutDir(resultD = resultD), cancer, "_", gene, "_", self, "_regulated_sites.txt"), quote = F, sep = "\t", row.names = F)
          
          ## make scatterplots
          phog <- loadPhosphoproteinNormalizedTumor(cancer = cancer)
          pho <- loadPhosphositeNormalizedTumor(cancer = cancer)
          pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho$Phosphosite), gene_vector = as.vector(pho$Gene))
          tab_gene <- tab_gene[order(tab_gene$FDR_pho_kin),]
          tab_gene$log10FDR <- -log10(tab_gene$FDR_pho_kin)
          tab_gene$circle_pt <- 5*tab_gene$log10FDR
          tab_gene$height_pt <- 10*tab_gene$log10FDR
          
          for (enzyme in unique(tab_gene$GENE)) {
            tab_en <- tab_gene[tab_gene$GENE == enzyme, ]
            tab_en$site <- paste0(tab_en$SUB_GENE, "_", tab_en$SUB_MOD_RSD)
            phog_en <- phog[phog$Gene == enzyme,]
            phog_en.m <- melt(phog_en)
            colnames(phog_en.m) <- c("GENE", "sampID", "phog_en")
            
            pho_sub <- cbind(pho, pho_head)[pho_head$SUBSTRATE == gene & pho_head$SUB_MOD_RSD %in% tab_gene$SUB_MOD_RSD[tab_gene$GENE == enzyme],]
            pho_sub.m <- melt(pho_sub)
            colnames(pho_sub.m) <- c("SUB_GENE", "Phosphosite", "transcript", "SUB_MOD_RSD", "substrate" , "sampID", "pho_sub")
            
            ## merge phosphorylation data together
            tab4plot <- merge(phog_en.m[, c("GENE", "sampID", "phog_en")], pho_sub.m[, c("SUB_GENE", "SUB_MOD_RSD", "sampID", "pho_sub")], by = c("sampID"), all = T)
            tab4plot$facet_y <- paste0(tab4plot$SUB_GENE, "_", tab4plot$SUB_MOD_RSD)
            tab4plot$subtype <- sampID2pam50(sampleID_vector = as.vector(tab4plot$sampID), pam50_map = pam50_map, sample_map = loadSampMap())
            tab4plot$facet_y <- factor(tab4plot$facet_y, levels = as.vector(tab_en$site))
            
            p <- ggplot()
            p <- p + geom_point(data = tab4plot, mapping = aes(x = phog_en, y = pho_sub, color = subtype), alpha = 0.8, shape = 16, stroke = 0)
            p <- p + facet_grid(facet_y~., scales = "fixed", space = "fixed")
            p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
            p <- p + scale_color_manual(values = brewer.pal(n = 6, name = "Set1"))
            p <- p + xlab(paste0(enzyme, " phosphorylation level")) + ylab(paste0(gene, " phosphosite level"))
            p <- p + theme(strip.text.y = element_text(size = 10, face = "bold"))
            p <- p + theme(axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"))
            p
            ggsave(filename = paste0(makeOutDir(resultD = resultD), enzyme, "_", gene, "_scatterplots.pdf"),
                   width = 4, height = pdf_sizes_per_yfacet[[as.character(length(unique(tab4plot$facet_y)))]])
          }
        }
      }
    }
    
    for (self in "cis") {
      for (cancer in "BRCA") {
        for (gene in "AKT1") {
          tab_gene <- sup_cans_tab_en_reg
          tab_gene <- tab_gene[tab_gene$SUB_GENE == gene & tab_gene$SELF == self & tab_gene$Cancer == cancer,]
          tab_gene <- merge(tab_gene, 
                            pho_head_cancers_sim[pho_head_cancers_sim$SUBSTRATE == gene & pho_head_cancers_sim$SUB_MOD_RSD %in% tab_gene$SUB_MOD_RSD, 
                                                 c("SUBSTRATE", "SUB_MOD_RSD", "seq_region_name", "start_hg19")], 
                            by.x = c("SUB_GENE", "SUB_MOD_RSD"), by.y = c("SUBSTRATE", "SUB_MOD_RSD"),
                            all.x = T)
          refseqs_mapTab = getBM(attributes = c("refseq_mrna", "refseq_peptide", "hgnc_symbol"),
                                 filters = c("refseq_peptide"),
                                 values = str_split_fixed(string = as.vector(tab_gene$transcript), pattern = "\\.", 2)[,1], 
                                 mart = ensembl, uniqueRows=T)
          
          
          tab4proteinpaint <- data.frame(gene = tab_gene$SUB_GENE, 
                                         refseq = refseqs_mapTab$refseq_mrna, 
                                         chromosome = tab_gene$seq_region_name, 
                                         start = tab_gene$start_hg19, 
                                         aachange = tab_gene$SUB_MOD_RSD, class = "silent")
          tab4proteinpaint <- tab4proteinpaint[!is.na(tab4proteinpaint$start),]
          tab4proteinpaint <- unique(tab4proteinpaint)
          write.table(x = tab4proteinpaint, file = paste0(makeOutDir(resultD = resultD), cancer, "_", gene, "_", self, "_regulated_sites.txt"), quote = F, sep = "\t", row.names = F)
          
          ## make scatterplots
          pro <- loadProteinNormalizedTumor(cancer = cancer)
          pho <- loadPhosphositeNormalizedTumor(cancer = cancer)
          pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho$Phosphosite), gene_vector = as.vector(pho$Gene))
          tab_gene <- tab_gene[order(tab_gene$FDR_pro_kin),]
          tab_gene$log10FDR <- -log10(tab_gene$FDR_pro_kin)
          tab_gene$circle_pt <- 3*tab_gene$log10FDR
          tab_gene$height_pt <- 10*tab_gene$log10FDR
          
          for (enzyme in unique(tab_gene$GENE)) {
            tab_en <- tab_gene[tab_gene$GENE == enzyme, ]
            tab_en$site <- paste0(tab_en$SUB_GENE, "_", tab_en$SUB_MOD_RSD)
            pro_en <- pro[pro$Gene == enzyme,]
            pro_en.m <- melt(pro_en)
            colnames(pro_en.m) <- c("GENE", "sampID", "pro_en")
            
            pho_sub <- cbind(pho, pho_head)[pho_head$SUBSTRATE == gene & pho_head$SUB_MOD_RSD %in% tab_gene$SUB_MOD_RSD[tab_gene$GENE == enzyme],]
            pho_sub.m <- melt(pho_sub)
            colnames(pho_sub.m) <- c("SUB_GENE", "Phosphosite", "transcript", "SUB_MOD_RSD", "substrate" , "sampID", "pho_sub")
            
            ## merge phosphorylation data together
            tab4plot <- merge(pro_en.m[, c("GENE", "sampID", "pro_en")], pho_sub.m[, c("SUB_GENE", "SUB_MOD_RSD", "sampID", "pho_sub")], by = c("sampID"), all = T)
            tab4plot$facet_y <- paste0(tab4plot$SUB_GENE, "_", tab4plot$SUB_MOD_RSD)
            tab4plot$subtype <- sampID2pam50(sampleID_vector = as.vector(tab4plot$sampID), pam50_map = pam50_map, sample_map = loadSampMap())
            tab4plot$facet_y <- factor(tab4plot$facet_y, levels = as.vector(tab_en$site))
            if (length(unique(tab4plot$facet_y)) > 5) {
              for (rsd in unique(tab4plot$facet_y)) {
                tab4plot_facet_y <- tab4plot[tab4plot$facet_y == rsd,]
                p <- ggplot()
                p <- p + geom_point(data = tab4plot_facet_y, mapping = aes(x = pro_en, y = pho_sub, color = subtype), alpha = 0.8, shape = 16, stroke = 0)
                p <- p + facet_grid(facet_y~., scales = "fixed", space = "fixed")
                p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
                p <- p + scale_color_manual(values = brewer.pal(n = 6, name = "Set1"))
                p <- p + xlab(paste0(enzyme, " protein level")) + ylab(paste0(gene, " phosphosite level"))
                p <- p + theme(strip.text.y = element_text(size = 10, face = "bold"))
                p <- p + theme(axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"))
                p
                ggsave(filename = paste0(makeOutDir(resultD = resultD), enzyme, "_", gene, "_", rsd, "_scatterplots.pdf"),
                       width = 4, height = pdf_sizes_per_yfacet[[as.character(length(unique(tab4plot_facet_y$facet_y)))]])
              }
            } else {
              p <- ggplot()
              p <- p + geom_point(data = tab4plot, mapping = aes(x = pro_en, y = pho_sub, color = subtype), alpha = 0.8, shape = 16, stroke = 0)
              p <- p + facet_grid(facet_y~., scales = "fixed", space = "fixed")
              p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
              p <- p + scale_color_manual(values = brewer.pal(n = 6, name = "Set1"))
              p <- p + xlab(paste0(enzyme, " protein level")) + ylab(paste0(gene, " phosphosite level"))
              p <- p + theme(strip.text.y = element_text(size = 10, face = "bold"))
              p <- p + theme(axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"))
              p
              ggsave(filename = paste0(makeOutDir(resultD = resultD), enzyme, "_", gene, "_scatterplots.pdf"),
                     width = 4, height = pdf_sizes_per_yfacet[[as.character(length(unique(tab4plot$facet_y)))]])
            }
          }
        }
      }
    }
    
  }
}