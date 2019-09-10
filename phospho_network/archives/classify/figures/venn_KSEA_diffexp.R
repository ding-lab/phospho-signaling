# Yige Wu @ WashU 2018 Apr
# venn diagram showing the agreement in the kinase activity inferred by substrate and by its own global phosphorylation


# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(ggplot2)
library(gplots)


# inputs ------------------------------------------------------------------

# generate venn diagrams ---------------------------------------------------
for (m.cutoff in 2:2 ) {
  for (NetworKIN.cutoff in 5:5) {
    ptms_site_pairs_sup <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv")
    ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$networkin_score >= NetworKIN.cutoff,]
    for (sig in c(0.2) ) {
      for (cancer in c("BRCA", "OV", "CO", "UCEC") ) {
        ## input kinase activity inferred by KSEA
        df <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA_genomic_matched/OmniPath_NetworKIN_", "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "/", cancer, "/KSEA Kinase Scores.csv"),
                    data.table = F)
        df <- df[df$p.value < 0.05 & df$m >= m.cutoff,]
        df$direction <- ifelse(df$z.score > 0, "up", "down")
        
        ## input differential expression for global phosphorylation
        pho <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", cancer, "_collapsed_PHO_diffexp_paired_FDR", sig, "_wilcoxon.txt"),
                     data.table = F)
        pho_kinase <- pho[pho$Gene %in% ptms_site_pairs_sup$GENE,]
        
        up_ksea <- as.vector(df$Kinase.Gene[df$direction == "up"])
        down_ksea <- as.vector(df$Kinase.Gene[df$direction == "down"])
        up_diff <- as.vector(pho_kinase$Gene[pho_kinase$diffexp_type == "up"])
        down_diff <- as.vector(pho_kinase$Gene[pho_kinase$diffexp_type == "down"])
        
        fn = paste0(makeOutDir(), "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "_", cancer, "_KSEA_diffexp_overlap.pdf")
        pdf(file = fn, width = 5, height = 5)
        venn(list(up_ksea=up_ksea,down_ksea=down_ksea,up_diff=up_diff,down_diff=down_diff))
        dev.off()
        
        sink(paste0(makeOutDir(), "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "_", cancer, "_KSEA_diffexp_overlap.txt"), append=FALSE, split=FALSE)
        for (ksea_dir in c("up", "down")) {
          for (diff_dir in c("up", "down")) {
            cat(paste0("overlap between KSEA ", ksea_dir, " kinases and differentially ", diff_dir, " phosphorylated kinases:\n"))
            cat(intersect(as.vector(df$Kinase.Gene[df$direction == ksea_dir]), as.vector(pho_kinase$Gene[pho_kinase$diffexp_type == diff_dir])))
            cat("\n")
          }
        }
        sink()
        closeAllConnections()
      }
    }
  }
}