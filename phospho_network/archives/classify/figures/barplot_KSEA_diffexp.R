# Yige Wu @ WashU 2018 Mar
# barplot displaying KSEA score for significantly enriched kinase/phosphotase

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

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
    # for (cancer in c("BRCA")) {
      df <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA/OmniPath_NetworKIN_", 
                                 "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "/", cancer, "/KSEA Kinase Scores.csv"),
                  data.table = F)
      df <- df[df$p.value < 0.05 & df$m >= m.cutoff,]
      df <- merge(df, op_enzyme, by.x = c("Kinase.Gene"), by.y = c("GENE"), all.x = T)
      df$enzyme_type[is.na(df$enzyme_type)] <- "kinase"
      df$direction <- ifelse(df$z.score > 0, "up", "down")
      ksea_diffexp <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/integrate_KSEA_diffexp/", 
                                           cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
      ksea_diffexp_direction <- data.frame(table(ksea_diffexp[, c("GENE", "substrate_direction")]))
      tmp <- merge(df, ksea_diffexp_direction, by.x = c("Kinase.Gene"), by.y = c("GENE"), all.x = T)
      tmp1 <- tmp[!is.na(tmp$substrate_direction),]
      tmp2 <- tmp[is.na(tmp$substrate_direction),]
      tmp3 <- tmp2; tmp3$substrate_direction <- "up"; tmp3$Freq <- 0
      tmp4 <- tmp3; tmp4$substrate_direction <- "down"
      df <- rbind(tmp1, tmp3, tmp4)
      df$Freq2p <- df$Freq; df$Freq2p[df$substrate_direction == "down"] <- (-df$Freq2p[df$substrate_direction == "down"])
      df$Kinase.Gene <- factor(df$Kinase.Gene, levels = as.vector(df$Kinase.Gene)[order(as.vector(df$z.score))])
      
      p <- ggplot()
      p <- p + geom_bar(data=df, aes(y = Freq2p, x = Kinase.Gene, fill = substrate_direction ), stat="identity")
      p <- p + geom_text(data=df, aes(y = ifelse(z.score > 0, -10, 10), x = Kinase.Gene, label = Kinase.Gene, color = direction ), size = 2)
      p <- p + facet_grid(.~enzyme_type, space = "free_y",scales = "free_y", drop = TRUE)#, space = "free", scales = "free")
      p <- p + scale_color_manual(values = c("up" = set1[1],"down" = set1[2]))
      p <- p + scale_fill_manual(values = c("up" = set1[1],"down" = set1[2]))
      p <- p + theme_bw() + theme_nogrid()
      p <- p + coord_flip()
      p <- p + xlab("enzyme")+ylab("number of differentially phosphorylated phosphosites")
      p <- p + theme(axis.title=element_text(size=10))
      p <- p + theme(axis.text.x = element_text(colour="black", size=10),
                     axis.text.y = element_blank())
      p <- p + ylim(min(-15, min(df$Freq2p)), max(15, max(df$Freq2p)))
      p
      fn = paste0(makeOutDir(), "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "_", cancer, "_KSEA_diffexp_substrates.pdf")
      ggsave(file=fn, height=10, width=6)
    }
  }
}