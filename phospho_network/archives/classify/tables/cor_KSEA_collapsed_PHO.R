# Yige Wu @ WashU 2018 Mar
# compare KSEA output with global phosphoprotein fold change between tumor and normal
# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/cptac2p_analysis_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')

coef_names <- c("correlation coefficient", "rho")
names(coef_names) <- c("pearson", "spearman")

for (method in c("pearson", "spearman")) {
  for (known_db in c("OmniPath", "PSP")) {
    for (addNetworKIN in c("", paste0("_NetworKIN", 2:5, "cutoff"))) {
      for (cancer in c("CO", "OV", "BRCA")) {
        phog_fc <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/generate_table4KSEA/", cancer, "_collapsed_PHO_FC_for_KSEA.csv"), data.table = F)
        ksea_score <- read.csv(file = paste0("cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA/", known_db, addNetworKIN, "/", cancer, "/KSEA Kinase Scores.csv"))
        df <- merge(ksea_score, phog_fc, by.x = c("Kinase.Gene"), by.y = c("Gene"), all.x = T)
        colnames(df)[ncol(df)] <- "Kinase.FC"
        df$Kinase.log2FC <- log2(df$Kinase.FC)
        
        ## do pearson correlation
        cor_stat <- cor.test(x = as.numeric(as.vector(df$z.score)), y = as.numeric(as.vector(df$Kinase.log2FC)), method = method)
        
        p = ggplot(df, aes(x=z.score, y=Kinase.log2FC))
        p = p + geom_point(alpha=0.5, stroke = 0)
        p = p + labs(x = paste0("PRM protein abundance"), 
                     y = paste0("global protein abundance"))
        p = p + annotate("text", x = quantile(df$z.score, probs = 0.95, na.rm = T), y = quantile(df$Kinase.log2FC, probs = 0.999, na.rm = T), label = paste0(method, "'s ", coef_names[method], " = ", signif(cor_stat$estimate[1], digits = 2)))
        p = p + annotate("text", x = quantile(df$z.score, probs = 0.95, na.rm = T), y = quantile(df$Kinase.log2FC, probs = 0.999, na.rm = T)-(range(df$Kinase.log2FC, na.rm = T)[2]-range(df$Kinase.log2FC, na.rm = T)[1])/20, label = paste0("p.value = ", signif(cor_stat$p.value[1], digits = 2)))
        p = p + theme_nogrid()
        p = p + theme(axis.title = element_text(size=10), 
                      axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), 
                      axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
        p = p + theme(title = element_text(size = 8))
        p
        resultDnow <- makeOutDir()
        fn = paste0(resultDnow, known_db, addNetworKIN, "_", cancer, "_", method, "_correlation_btw_KSEA_and_kinase_collapsed_PHO.pdf")
        ggsave(file=fn, height=5, width=5, useDingbats=FALSE)
      }
    }
  }
  
}