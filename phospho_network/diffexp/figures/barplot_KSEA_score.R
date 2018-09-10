# Yige Wu @ WashU 2018 Mar
# barplot displaying KSEA score for significantly enriched kinase/phosphotase

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')

# test --------------------------------------------------------------------
## input omnipath to annotate kinase/phosphotase
op <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/compile_omnipath/omnipath_genename_annotated.csv")
op$enzyme_type <- NA
op$enzyme_type[op$modification == "phosphorylation"] <- "kinase"
op$enzyme_type[op$modification == "dephosphorylation"] <- "phosphotase"
op_enzyme <- unique(op[op$enzyme_type == "kinase" | op$enzyme_type == "phosphotase", c("GENE", "enzyme_type")])

for (m.cutoff in 2:2) {
  for (NetworKIN.cutoff in 5:5) {
    for (cancer in c("BRCA", "OV", "CO")) {
      df <- fread(input = paste0("/Users/yigewu/Box\ Sync/", "cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA/OmniPath_NetworKIN_", 
                                 "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "/", cancer, "/KSEA Kinase Scores.csv"),
                  data.table = F)
      df <- df[df$p.value < 0.05 & df$m >= m.cutoff,]
      df <- merge(df, op_enzyme, by.x = c("Kinase.Gene"), by.y = c("GENE"), all.x = T)
      df$enzyme_type[is.na(df$enzyme_type)] <- "kinase"
      df$direction <- ifelse(df$z.score > 0, "up", "down")
      df$Kinase.Gene <- factor(df$Kinase.Gene, levels = as.vector(df$Kinase.Gene)[order(as.vector(df$z.score))])
      
      p <- ggplot()
      p <- p + geom_bar(data=df, aes(y = z.score, x = Kinase.Gene, fill = direction ), stat="identity")
      p <- p + geom_text(data=df, aes(y = ifelse(z.score > 0, -1.5, 1.5), x = Kinase.Gene, label = Kinase.Gene, color = direction), size = 2)
      p <- p + facet_grid(.~enzyme_type, space = "free_y",scales = "free_y", drop = TRUE)#, space = "free", scales = "free")
      p <- p + scale_color_manual(values = c("up" = set1[1],"down" = set1[2]))
      p <- p + scale_fill_manual(values = c("up" = set1[1],"down" = set1[2]))
      p <- p + theme_bw() + theme_nogrid()
      p <- p + coord_flip()
      p <- p + xlab("enzyme")+ylab("KSEA z.score")
      p <- p + theme(axis.title=element_text(size=10))
      p <- p + theme(axis.text.x = element_text(colour="black", size=10), 
                     axis.text.y = element_blank())
      p
      fn = paste0(makeOutDir(), "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "_", cancer, "_KSEA_z.score.pdf")
      ggsave(file=fn, height=10, width=6)
    }
  }
}


