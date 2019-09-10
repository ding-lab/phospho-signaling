# Yige Wu @ WashU 2018 Apr
# barplot displaying KSEA score for significantly enriched kinase/phosphotase

# source ------------------------------------------------------------------
# source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(ggplot2)
# test --------------------------------------------------------------------
## input omnipath to annotate kinase/phosphotase
op <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/compile_omnipath/omnipath_genename_annotated.csv")
op$enzyme_type <- NA
op$enzyme_type[op$modification == "phosphorylation"] <- "kinase"
op$enzyme_type[op$modification == "dephosphorylation"] <- "phosphotase"
op_kinase <- unique(op[op$enzyme_type == "kinase" & !is.na(op$enzyme_type), c("GENE", "enzyme_type")])
op_phosphotase <- unique(op[op$enzyme_type == "phosphotase" & !is.na(op$enzyme_type), c("GENE", "enzyme_type")])
op_enzyme <- rbind(op_phosphotase, op_kinase[!(op_kinase$GENE %in% op_phosphotase$GENE),])


for (m.cutoff in 2:2 ) {
  for (NetworKIN.cutoff in 5:5) {
    for (sig in c(0.2) ) {
      for (cancer in c("BRCA", "OV", "CO", "UCEC") ) {
        df <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA_genomic_matched/OmniPath_NetworKIN_", "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "/", cancer, "/KSEA Kinase Scores.csv"),
                    data.table = F)
        df <- df[df$p.value < 0.05 & df$m >= m.cutoff,]
        df <- merge(df, op_enzyme, by.x = c("Kinase.Gene"), by.y = c("GENE"), all.x = T)
        df$enzyme_type[is.na(df$enzyme_type)] <- "kinase"
        df$direction <- ifelse(df$z.score > 0, "up", "down")
        
        ## input differential expression
        pho <- fread(paste0(resultD, "diffexp/tables/differential_expression_paired/", cancer, "_PHO_diffexp_paired_FDR", sig, "_wilcoxon.txt"),
                     data.table = F)
        pho <- cbind(pho, formatPhosphosite(phosphosite_vector = pho$Phosphosite, pho$Gene))
        pho <- pho[pho$Gene %in% df$Kinase.Gene,]
        
        ## add in the fold change for all substrates
        fc <- fread(input = paste0("/Users/yigewu/Box\ Sync/", "cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/generate_table4KSEA/", cancer, "_phosphposite_FC_for_KSEA.csv"),
                    data.table = F)
        fc <- fc[!is.na(fc$FC) & fc$Gene %in% df$Kinase.Gene,]
        fc$direction <- ifelse(fc$FC > 1, "up", "down")
        fc <- merge(fc, pho, by.x = c("Gene", "Residue.Both"), by.y = c("Gene", "SUB_MOD_RSD"), all.x = T)
        fc$diffexp <- ifelse(!is.na(fc$`q-value(%)`), paste0("FDR<", sig), paste0("FDR>=", sig))
        fc$direction_diffexp <- paste0(fc$direction, " (", fc$diffexp, ")")
        fc_kinase <- data.frame(table(fc[, c("Gene", "direction_diffexp")]))
        
        df <- merge(df, fc_kinase, by.x = c("Kinase.Gene"), by.y = c("Gene"), all.x = T)
        tmp <- as.vector(df$Freq)
        tmp[is.na(tmp)] <- 0
        df$Freq <- tmp
        tmp <- as.vector(df$direction_diffexp)
        tmp[is.na(tmp)] <- "down (FDR<0.2)"
        df$direction_diffexp <- tmp
        tmp <- as.vector(df$Freq)
        tmp[grepl(x = df$direction_diffexp, pattern = "down")] <- -(tmp[grepl(x = df$direction_diffexp, pattern = "down")])
        df$Freq2p <- tmp
        df$substrate_direction <- ifelse(grepl(pattern = "up", x = df$direction_diffexp), "up", "down")
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
        df$direction_diffexp <- as.vector(df$direction_diffexp)
        tmp <- unique(df[,c("Kinase.Gene", "z.score")])
        df$Kinase.Gene <- factor(df$Kinase.Gene, levels = as.vector(tmp$Kinase.Gene)[order(as.vector(tmp$z.score))])
        
        
        p <- ggplot()
        p <- p + geom_bar(data=df, aes(y = Freq2p, x = Kinase.Gene, fill = direction_diffexp ), stat="identity")
        p <- p + geom_text(data=df, aes(y = x, x = Kinase.Gene, label = Kinase.Gene), size = 2, color = "black" )
        p <- p + scale_fill_manual(values = c("down (FDR<0.2)" = "#1F78B4", "up (FDR<0.2)" = "#E31A1C", 
                                              "down (FDR>=0.2)" = "#A6CEE3", "up (FDR>=0.2)" = "#FB9A99"))
        p <- p + theme_bw() + theme_nogrid()
        p <- p + facet_grid(~enzyme_type, scales = "free")
        p <- p + coord_flip()
        p <- p + xlab("enzyme")+ylab("number of differentially phosphorylated phosphosites on the kinase")
        p <- p + theme(axis.title=element_text(size=10))
        p <- p + theme(axis.text.x = element_text(colour="black", size=10),
                       axis.text.y = element_blank())
        p
        fn = paste0(makeOutDir(), "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "_", cancer, "_KSEA_diffexp_substrates.pdf")
        ggsave(file=fn, height=10, width=6)
      }
    }
  }
}





