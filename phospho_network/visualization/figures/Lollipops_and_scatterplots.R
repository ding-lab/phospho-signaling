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
reg_nonNA2test <- c(20)
top_kinase2show <- 10
top_substrate2show <- 2
pdf_sizes <- list(trans = list(BRCA = list(width = 7, height = 6),
                               OV = list(width = 4, height = 3),
                               CO = list(width = 5, height = 4)))
pdf_sizes_per_yfacet <- list("1" = 2.5, "2" = 4, "3" = 5, "4" = 6, "5" = 7, "16" = 20, "20" = 25)
clinical <- loadSampMap()
nrow(clinical)
fdr_thres2process <- c(0.05)
cptac_phase2process <- "cptac2p_3can"
# inputs -------------------------------------------------------------------
## input driver genes
# driver_genes <- read.delim(file = "./TCGA_data/reference_files/Consensus.Genelist.full.txt")
# pho_head_cancers_sim <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid.txt"), data.table = F)
pho_head_cancers_sim <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid_hg38_hg19.txt"), data.table = F)

driver_pancan <- loadGeneList(gene_type = "driver", cancer = "PANCAN", is.soft.limit = "")

ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org") ## 38
filters = listFilters(ensembl)

pam50_map <- loadPAM50Map()

## input phosphorylation site with uniprots
op <- fread(input = "cptac2p/analysis_results/preprocess_files/tables/convert_omnipath2genomic_coord/op_w.phosphosite.availability.txt", data.table = F)
op$SUB_MOD_RSD <- paste0(op$residue_type, op$residue_offset)

## input CPTAC ID mapping
map_protein_id_map <- fread(input = paste0(cptac_sharedD, "UniProt.human.20171228.RISnrNF.259contams.fasta.crossRefs.RefSeq.20171003_Human_ucsc_hg38_cpdb_mito_259contamsnr.fasta.categories.tsv"),
                            data.table = F)
## annotate kinase substrate regulation
for (enzyme_type in c("kinase")) {
  subdir1 <- paste0(makeOutDir(resultD = resultD), enzyme_type, "/")
  dir.create(subdir1)
  for (reg_nonNA in reg_nonNA2test) {
    # reg_nonNA <- 15
    subdir2 <- paste0(subdir1, "reg_nonNA", reg_nonNA, "/")
    dir.create(subdir2)
    
    sup_cans_tab_en <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                            enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor",
                                            "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    for (fdr_thres in fdr_thres2process) {
      # fdr_thres <- 0.3
      subdir3 <- paste0(subdir2, "fdr_thres", fdr_thres, "/")
      dir.create(subdir3)
      
      sup_cans_tab_en <- markSigSiteCan(regression = sup_cans_tab_en, sig_thres = fdr_thres, enzyme_type = enzyme_type)
      sup_cans_tab_en$regulated <- (sup_cans_tab_en$fdr_sig & sup_cans_tab_en$coef_sig)
      # sup_cans_tab_en_reg <- sup_cans_tab_en[sup_cans_tab_en$regulated,]
      
      for (self in "cis") {
        # tab_self <- sup_cans_tab_en_reg
        tab_self <- sup_cans_tab_en
        tab_self <- tab_self[tab_self$SELF == self, ]
        if (enzyme_type == "kinase") {
          tab_self_path <- tab_self[tab_self$GENE %in% names(unlist(map2TCGApathwaways(gene_list = as.vector(tab_self$SUB_GENE), pathway_list = tcga_pathways))),]
        }
        tab_self_path <- tab_self

        for (gene in c("ERBB2", "AKT1", "GSK3B", "BRAF", "EGFR", "IGF1R", "MET")) {
        # for (gene in c("ERBB2")) {
          subdir4 <- paste0(subdir3, gene, "/")
          dir.create(subdir4)
          
          tab_gene <- tab_self_path
          tab_gene <- tab_gene[tab_gene$SUB_GENE == gene,]
          tab_gene <- merge(tab_gene, map_protein_id_map[, c( "accession_number", "geneSymbol")], 
                            by.x = c("SUB_GENE"), by.y = c("geneSymbol"),
                            all.x = T)
          tab_gene$accession_number <- str_split_fixed(string = tab_gene$accession_number, pattern = "-", 2)[,1]
          tab_gene <- unique(tab_gene)
          
          tab_gene <- tab_gene[order(tab_gene$FDR_pro_kin),]
          tab_gene$log10FDR <- -log10(tab_gene$FDR_pro_kin)
          tab_gene$circle_pt <- signif(2^tab_gene$log10FDR, digits = 0)
          tab_gene$height_pt <- signif(10*tab_gene$log10FDR, digits = 0)
          
          ## make commands for lollipops
          uniprotID <- as.character(unique(tab_gene$accession_number[!is.na(tab_gene$accession_number)]))
          print(uniprotID)
          uniprotID_sims <- str_split_fixed(uniprotID, pattern = "-", 2)[,1]
          for (uniprotID_sim in uniprotID_sims) {
            cmd <- paste0("cp -r ./cptac2p_analysis/dependencies/lollipops_1.4.0_mac64 ", subdir4, ";cd ", subdir4, "; ./lollipops_1.4.0_mac64/lollipops -legend -labels -U ", uniprotID_sim, 
                          " -o=", gene, "_", uniprotID_sim, ".png -dpi=800 -domain-labels=fit")
            for (rsd in unique(tab_gene$SUB_MOD_RSD)) {
              if (str_count(rsd, "[STY]") == 1) {
                aa <- substr(x = rsd, start = 1, stop = 1)
                cmd <- paste0(cmd, " ", rsd, color_aa[aa], "@", tab_gene$circle_pt[tab_gene$SUB_MOD_RSD == rsd])
              }
            }
            system(command = cmd)
          }
          
          for (uniprotID_sim in uniprotID_sims) {
            cmd <- paste0("cp -r ./cptac2p_analysis/dependencies/lollipops_1.4.0_mac64 ", subdir4, ";cd ", subdir4, "; ./lollipops_1.4.0_mac64/lollipops -legend -labels -U ", uniprotID_sim, 
                          " -o=", gene, "_", uniprotID_sim,"_FDR", fdr_thres, ".png -dpi=800 -domain-labels=fit")
            for (rsd in unique(tab_gene$SUB_MOD_RSD[tab_gene$regulated])) {
              if (str_count(rsd, "[STY]") == 1) {
                aa <- substr(x = rsd, start = 1, stop = 1)
                cmd <- paste0(cmd, " ", rsd, color_aa[aa], "@", tab_gene$circle_pt[tab_gene$SUB_MOD_RSD == rsd])
              }
            }
            system(command = cmd)
          }

          
          # for (cancer in unique(tab_gene$Cancer)) {
          #   ## make scatterplots
          #   if (self == "cis") {
          #     pro <- loadProteinNormalizedTumor(cancer = cancer)
          #     pho <- loadPhosphositeNormalizedTumor(cancer = cancer)
          #     pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho$Phosphosite), gene_vector = as.vector(pho$Gene))
          #     
          #     for (enzyme in unique(tab_gene$GENE)) {
          #       tab_en <- tab_gene[tab_gene$GENE == enzyme & tab_gene$Cancer == cancer, ]
          #       tab_en$site <- paste0(tab_en$SUB_GENE, "_", tab_en$SUB_MOD_RSD)
          #       pro_en <- pro[pro$Gene == enzyme,]
          #       pro_en.m <- melt(pro_en)
          #       colnames(pro_en.m) <- c("GENE", "sampID", "pro_en")
          #       
          #       pho_sub <- cbind(pho, pho_head)[pho_head$SUBSTRATE == gene & pho_head$SUB_MOD_RSD %in% tab_gene$SUB_MOD_RSD[tab_gene$GENE == enzyme],]
          #       pho_sub.m <- melt(pho_sub)
          #       colnames(pho_sub.m) <- c("SUB_GENE", "Phosphosite", "transcript", "SUB_MOD_RSD", "substrate" , "sampID", "pho_sub")
          #       
          #       ## merge phosphorylation data together
          #       tab4plot <- merge(pro_en.m[, c("GENE", "sampID", "pro_en")], pho_sub.m[, c("SUB_GENE", "SUB_MOD_RSD", "sampID", "pho_sub")], by = c("sampID"), all = T)
          #       if (cancer == "BRCA") {
          #         tab4plot$subtype <- sampID2pam50(sampleID_vector = as.vector(tab4plot$sampID), pam50_map = pam50_map, sample_map = loadSampMap())
          #       } else if (cancer == "CO") {
          #         tab4plot$subtype <- sampID2MSI(sampleID_vector = as.vector(tab4plot$sampID), subtype_map = loadMSIMap(), sample_map = loadSampMap())
          #       } else {
          #         tab4plot$subtype <- ""
          #       }
          #       for (rsd in unique(tab_en$SUB_MOD_RSD)) {
          #         tab4plot_facet_y <- tab4plot[tab4plot$SUB_MOD_RSD == rsd,]
          #         tab4plot_facet_y$facet_y <- paste0(gene, "_", rsd)
          #         p <- ggplot()
          #         p <- p + geom_point(data = tab4plot_facet_y, mapping = aes(x = pro_en, y = pho_sub, color = subtype), alpha = 0.8, shape = 16, stroke = 0)
          #         p <- p + facet_grid(facet_y~., scales = "fixed", space = "fixed")
          #         p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
          #         p <- p + scale_color_manual(values = brewer.pal(n = 6, name = "Set1"))
          #         p <- p + xlab(paste0(enzyme, " protein level (", cancer, ")")) + ylab(paste0(gene, " phosphosite level"))
          #         p <- p + theme(strip.text.y = element_text(size = 10, face = "bold"))
          #         p <- p + theme(axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"))
          #         p
          #         ggsave(filename = paste0(subdir4, cancer, "_", enzyme, "_", gene, "_", rsd, "_scatterplots.pdf"),
          #                width = 4, height = pdf_sizes_per_yfacet[[as.character(length(unique(tab4plot_facet_y$facet_y)))]])
          #       }
          #     }
          #   }
          #   if (self == "trans") {
          #     phog <- loadPhosphoproteinNormalizedTumor(cancer = cancer)
          #     pho <- loadPhosphositeNormalizedTumor(cancer = cancer)
          #     pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho$Phosphosite), gene_vector = as.vector(pho$Gene))
          #     
          #     for (enzyme in unique(tab_gene$GENE[tab_gene$Cancer == cancer & tab_gene$regulated])) {
          #       tab_en <- tab_gene[tab_gene$GENE == enzyme & tab_gene$Cancer == cancer, ]
          #       tab_en$site <- paste0(tab_en$SUB_GENE, "_", tab_en$SUB_MOD_RSD)
          #       phog_en <- phog[phog$Gene == enzyme,]
          #       phog_en.m <- melt(phog_en)
          #       colnames(phog_en.m) <- c("GENE", "sampID", "pho_en")
          #       
          #       pho_sub <- cbind(pho, pho_head)[pho_head$SUBSTRATE == gene & pho_head$SUB_MOD_RSD %in% tab_gene$SUB_MOD_RSD[tab_gene$GENE == enzyme],]
          #       pho_sub.m <- melt(pho_sub)
          #       colnames(pho_sub.m) <- c("SUB_GENE", "Phosphosite", "transcript", "SUB_MOD_RSD", "substrate" , "sampID", "pho_sub")
          #       
          #       ## merge phosphorylation data together
          #       tab4plot <- merge(phog_en.m[, c("GENE", "sampID", "pho_en")], pho_sub.m[, c("SUB_GENE", "SUB_MOD_RSD", "sampID", "pho_sub")], by = c("sampID"), all = T)
          #       if (cancer == "BRCA") {
          #         tab4plot$subtype <- sampID2pam50(sampleID_vector = as.vector(tab4plot$sampID), pam50_map = pam50_map, sample_map = loadSampMap())
          #       } else if (cancer == "CO") {
          #         tab4plot$subtype <- sampID2MSI(sampleID_vector = as.vector(tab4plot$sampID), subtype_map = loadMSIMap(), sample_map = loadSampMap())
          #       } else {
          #         tab4plot$subtype <- ""
          #       }
          #       for (rsd in unique(tab_en$SUB_MOD_RSD[tab_en$regulated])) {
          #         tab4plot_facet_y <- tab4plot[tab4plot$SUB_MOD_RSD == rsd,]
          #         tab4plot_facet_y$facet_y <- paste0(gene, "_", rsd)
          #         p <- ggplot()
          #         p <- p + geom_point(data = tab4plot_facet_y, mapping = aes(x = pho_en, y = pho_sub, color = subtype), alpha = 0.8, shape = 16, stroke = 0)
          #         p <- p + facet_grid(facet_y~., scales = "fixed", space = "fixed")
          #         p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
          #         p <- p + scale_color_manual(values = brewer.pal(n = 6, name = "Set1"))
          #         p <- p + xlab(paste0(enzyme, " phosphorylation level (", cancer, ")")) + ylab(paste0(gene, " phosphosite level"))
          #         p <- p + theme(strip.text.y = element_text(size = 10, face = "bold"))
          #         p <- p + theme(axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"))
          #         p
          #         ggsave(filename = paste0(subdir4, cancer, "_", enzyme, "_", gene, "_", rsd, "_scatterplots.pdf"),
          #                width = 4, height = pdf_sizes_per_yfacet[[as.character(length(unique(tab4plot_facet_y$facet_y)))]])
          #       }
          #     }
          #   }
          # }
        }
      }
    }
  }
}

