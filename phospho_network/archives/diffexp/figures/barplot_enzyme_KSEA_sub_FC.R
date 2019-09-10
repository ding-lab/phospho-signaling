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
op_kinase <- unique(op[op$enzyme_type == "kinase" & !is.na(op$enzyme_type), c("GENE", "enzyme_type")])
op_phosphotase <- unique(op[op$enzyme_type == "phosphotase" & !is.na(op$enzyme_type), c("GENE", "enzyme_type")])
op_enzyme <- rbind(op_phosphotase, op_kinase[!(op_kinase$GENE %in% op_phosphotase$GENE),])


for (m.cutoff in 2:2) {
  for (NetworKIN.cutoff in 5:5) {
    for (cancer in c("BRCA", "OV", "CO", "UCEC")) {
    # for (cancer in c("BRCA")) {
      df <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA_genomic_matched/OmniPath_NetworKIN_", 
                                 "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "/", cancer, "/KSEA Kinase Scores.csv"),
                  data.table = F)
      df <- df[df$p.value < 0.05 & df$m >= m.cutoff,]
      df <- merge(df, op_enzyme, by.x = c("Kinase.Gene"), by.y = c("GENE"), all.x = T)
      df$enzyme_type[is.na(df$enzyme_type)] <- "kinase"
      df$direction <- ifelse(df$z.score > 0, "up", "down")
      ksea_diffexp <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/diffexp/tables/integrate_KSEA_diffexp_FC/", 
                                           cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
      ksea_diffexp_direction <- data.frame(table(ksea_diffexp[, c("GENE", "substrate_FC_diffexp")]))
      df <- merge(df, ksea_diffexp_direction, by.x = c("Kinase.Gene"), by.y = c("GENE"), all.x = T)
      df$Freq2p <- df$Freq; df$Freq2p[grepl(x = df$substrate_FC_diffexp, pattern = "down")] <- (-df$Freq2p[grepl(x = df$substrate_FC_diffexp, pattern = "down")])
      df$Freq2p <- as.numeric(df$Freq2p)
      df$substrate_direction <- ifelse(grepl(pattern = "up", x = df$substrate_FC_diffexp), "up", "down")
      x <- sapply(unique(df$Kinase.Gene), FUN = function(g, dat) {
        rows <- which(dat$Kinase.Gene == g)
        direction_g <- unique(dat$direction[rows])
        et <- unique(dat$enzyme_type[rows])
        if (direction_g == "up") {
          xs <- 0.5*(min(dat$Freq2p[dat$enzyme_type == et & dat$substrate_direction != direction_g])) + 0.3
        } else {
          xs <- 0.5*(max(dat$Freq2p[dat$enzyme_type == et & dat$substrate_direction != direction_g])) + 0.3
        }
        return(xs)
      }, dat = df)
      
      tmp <- data.frame(Kinase.Gene = unique(df$Kinase.Gene), x = x)
      df <- merge(df, tmp, all.x = T)
      df$Kinase.Gene <- factor(df$Kinase.Gene, levels = as.vector(df$Kinase.Gene)[order(as.vector(df$z.score))])

      p <- ggplot()
      p <- p + geom_bar(data=df, aes(y = Freq2p, x = Kinase.Gene, fill = substrate_FC_diffexp ), stat="identity")
      p <- p + geom_text(data=df, aes(y = x, x = Kinase.Gene, label = Kinase.Gene), size = 2, color = "black" )
      p <- p + scale_fill_manual(values = c("down (FDR<0.2)" = "#1F78B4", "up (FDR<0.2)" = "#E31A1C", "down (FDR>=0.2)" = "#A6CEE3", "up (FDR>=0.2)" = "#FB9A99"))
      p <- p + theme_bw() + theme_nogrid()
      p <- p + facet_grid(~enzyme_type, scales = "free")
      p <- p + coord_flip()
      p <- p + xlab("enzyme")+ylab("number of differentially phosphorylated phosphosites")
      p <- p + theme(axis.title=element_text(size=10))
      p <- p + theme(axis.text.x = element_text(colour="black", size=10),
                     axis.text.y = element_blank())
      p
      fn = paste0(makeOutDir(), "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "_", cancer, "_KSEA_diffexp_substrates.pdf")
      ggsave(file=fn, height=10, width=6)
    }
  }
}
