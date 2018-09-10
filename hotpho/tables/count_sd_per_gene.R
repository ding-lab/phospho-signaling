# Yige Wu @ WashU 2018 Aug
## generate table for plotting lolliplots for hotpho hybrid ptms
## number of observations and SD


# source ------------------------------------------------------------------
setwd("~/Box Sync")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")

# inputs hotpho results (hg38 based)------------------------------------------------------------------
hotpho <- read_delim("~/Box Sync/Ding_Lab/Projects_Current/hotpho_data/output/Data_201807_cc.p0.05.cluster_transcriptSynced_hybrid_SourceAnnotated.tsv", 
                     "\t", escape_double = FALSE, col_types = cols(Cluster = col_character()), trim_ws = TRUE)
hotpho <- data.frame(hotpho)

driver_table <- read_excel("./Ding_Lab/Projects_Current/TCGA_data/gene_lists/mmc1.xlsx", 
                           sheet = "Table S1", skip = 3)
driver_table <- data.frame(driver_table)

# take hotpho clusters with ptms annotated to CPTAC dataset and split per cohort---------------
hotpho_cptac <- hotpho[!is.na(hotpho$Cohort) & hotpho$Cohort != "Not_Found",]
table(hotpho_cptac$Cohort)
cohort_list <- sapply(hotpho_cptac$Cohort, FUN = function(c) strsplit(x = c, split = "\\|")[[1]], simplify = T)
cohort <- unlist(cohort_list, use.names = F)
cohort_num <- sapply(cohort_list, FUN = function(c) length(c))
hotpho_cptac_m <- hotpho_cptac[rep(x = 1:nrow(hotpho_cptac), cohort_num),]
hotpho_cptac_m$Cohort <- cohort
## remove the retrospective ones
hotpho_cptac_m <- hotpho_cptac_m[!grepl(pattern = "Retrospective", x = hotpho_cptac_m$Cohort) & !grepl(pattern = "CPTAC3", x = hotpho_cptac_m$Cohort),]
hotpho_cptac_m$Cohort <- as.vector(str_split_fixed(string = hotpho_cptac_m$Cohort, pattern = "-", n = 2)[, 2])
hotpho_cptac_m$rsd <- str_split_fixed(string = hotpho_cptac_m$Mutation_Gene, pattern = "\\.", 2)[, 2]
write.table(x = hotpho_cptac_m, file = paste0(makeOutDir(resultD = resultD), "hotpho_ptms_cptac2_prospective.txt"), quote = F, sep = '\t', row.names = F)

# set variables -----------------------------------------------------------
## assign colors to each cohorts
color_cohorts <- brewer.pal(length(unique(hotpho_cptac_m$Cohort)), "Set1")
names(color_cohorts) <- unique(hotpho_cptac_m$Cohort)
judeClasses = c('M','E','F','N','S','D','I','P','L','Intron','Utr3','Utr5','X','noncoding','snv','mnv','insertion','deletion')
cohorts2process <- c("BRCA", "OV", "CO", "CCRCC", "UCEC", "LUAD")
# genes2test <- c("BRAF", "PIK3R1", "CTNNB1")
genes2test <- c("BRAF", "PIK3R1", "CTNNB1", "AKT1", "CDK12", "TP53")
# genes2test <- unique(hotpho_cptac_m$Gene_Drug)
hotpho_cptac_m_genes <- NULL
for (gene in genes2test) {
  hotpho_cptac_m_g <- hotpho_cptac_m[hotpho_cptac_m$Gene_Drug == gene,]

  hotpho_cptac_m_g_new <- NULL
  for (cohort in unique(hotpho_cptac_m_g$Cohort)) {
    # cohort <- "BRCA"
    if (cohort != "LUAD") {
      pho <- loadPhosphositeNormalizedTumor(cancer = cohort)
      pho <- pho[pho$Gene == gene,]
      pho$rsd <- toupper(str_split_fixed(string = pho$Phosphosite, pattern = ":", 2)[,2])
      samples <- colnames(pho); samples <- samples[!(samples %in% c("Gene", "Phosphosite"))]
      for (rsd in unique(hotpho_cptac_m_g$rsd[hotpho_cptac_m_g$Cohort == cohort])) {
        # rsd <- "S681"
        pho_rsd <- pho[grepl(pattern = rsd, x = pho$Phosphosite, ignore.case = T), samples]
        if (nrow(pho_rsd) == 0) {
          stop(paste0("manual check ", rsd))
        } else if (nrow(pho_rsd) == 1){
          row2edit <- (hotpho_cptac_m_g$Cohort == cohort & hotpho_cptac_m_g$rsd == rsd)
          hotpho_cptac_m_g2add <- hotpho_cptac_m_g[row2edit,]
          hotpho_cptac_m_g2add$count_sample <- rowSums(x = !is.na(pho_rsd), na.rm = T)
          hotpho_cptac_m_g2add$sd_sample <- sd(x = as.vector(pho_rsd), na.rm = T)
          hotpho_cptac_m_g_new <- rbind(hotpho_cptac_m_g_new, hotpho_cptac_m_g2add)
          
        } else if (nrow(pho_rsd) > 1) {
          row2edit <- (hotpho_cptac_m_g$Cohort == cohort & hotpho_cptac_m_g$rsd == rsd)
          hotpho_cptac_m_g2add <- hotpho_cptac_m_g[row2edit,]
          for (rsd2 in pho$rsd[grepl(pattern = rsd, x = pho$Phosphosite, ignore.case = T)]) {
            hotpho_cptac_m_g2add2 <- hotpho_cptac_m_g2add
            hotpho_cptac_m_g2add2$rsd <- rsd2
            
            pho_rsd2 <- pho[pho$rsd == rsd2, samples]
            pho_rsd2 <- pho_rsd2[1,]
            hotpho_cptac_m_g2add2$count_sample <- rowSums(x = !is.na(pho_rsd2), na.rm = T)
            hotpho_cptac_m_g2add2$sd_sample <- sd(x = as.vector(pho_rsd2), na.rm = T)
            
            hotpho_cptac_m_g_new <- rbind(hotpho_cptac_m_g_new, hotpho_cptac_m_g2add2)
          }
        }
       }
    }
  }
  hotpho_cptac_m_g_new$judeClass <- sapply(X = hotpho_cptac_m_g_new$Cohort, FUN = function(c, jcs, cs) jcs[which(cs == c)], jcs = judeClasses, cs = cohorts2process)
  hotpho_cptac_m_g_new$text4proteinpaint <- paste0(hotpho_cptac_m_g_new$Mutation_Gene, "(", hotpho_cptac_m_g_new$Cohort, ")")
  hotpho4jude <- hotpho_cptac_m_g_new[rep(1:nrow(hotpho_cptac_m_g_new), as.vector(hotpho_cptac_m_g_new$count_sample)), c("text4proteinpaint", "Position", "judeClass")]
  write.table(x = hotpho4jude, file = paste0(makeOutDir(resultD = resultD), "hotpho4jude_", gene,  ".txt"), row.names = F, quote = F, sep = ";", col.names = F)
  hotpho_cptac_m_genes <- rbind(hotpho_cptac_m_genes, hotpho_cptac_m_g_new)
}


# plot bubble plots for SDs and number of samples detected ----------------
hotpho_cptac_m_genes <- data.frame(hotpho_cptac_m_genes)
library(ggrepel)
tab4plot <- hotpho_cptac_m_genes
tab4plot <- tab4plot[tab4plot$sd_sample > 0 & !is.na(tab4plot$sd_sample),]
tab4plot_sorted <- NULL
for (gene in unique(tab4plot$Gene_Drug)) {
  tab_gene <- tab4plot[tab4plot$Gene_Drug == gene,]
  tab_gene <- tab_gene[order(tab_gene$Position),]
  positions <- unique(as.vector(tab_gene$Position))
  for (pos in positions) {
    tab_pos <- tab_gene[tab_gene$Position == pos,]
    tab4plot_sorted <- rbind(tab4plot_sorted, tab_pos)
  }
}
tab4plot_sorted$Gene_Drug <- factor(tab4plot_sorted$Gene_Drug, levels = unique(tab4plot_sorted$Gene_Drug))
tab4plot_sorted$Mutation_Gene <- factor(tab4plot_sorted$Mutation_Gene, levels = rev(unique(tab4plot_sorted$Mutation_Gene)))

p <- ggplot(tab4plot_sorted) 
p <- p + geom_point(aes(x = Mutation_Gene, y = Cohort, color = Cohort, size = sd_sample), 
                    shape = 16, stroke = 0, alpha = 0.6)
p <- p + facet_grid(Gene_Drug~., scales = "free", space = "free")
p <- p + theme_bw()
p <- p + theme(axis.text.y = element_text(colour="black", size=10))
p <- p + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 0.5, vjust = 0.5))
p <- p + theme(strip.background.x = element_rect(color = "black"))
p <- p + theme(panel.spacing.y = unit(0, "lines"), strip.text.y = element_text(size = 8, angle = 0, face = "bold"))
p <- p + coord_flip()
p
fn = paste0(makeOutDir(resultD = resultD), "sample_sd.pdf")
ggsave(file=fn, height=5, width=3, useDingbats=FALSE)





