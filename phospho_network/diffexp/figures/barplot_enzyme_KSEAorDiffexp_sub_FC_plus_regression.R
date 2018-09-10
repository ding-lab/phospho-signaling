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
    ## input enzyme-substrate table
    ptms_site_pairs_sup <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv")
    ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$networkin_score >= NetworKIN.cutoff,]
    for (cancer in c("BRCA", "OV", "CO", "UCEC")) {
      for (reg_sig in c(0.05)) {
        ## get KSEA scores integrated with differential phosphorylation results
        ksea_diffexp <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/diffexp/tables/integrate_enzyme_KSEAordiffexp_sub_FC_plus_regression/", 
                                             cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
        
        ksea_diffexp2 <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/diffexp/tables/integrate_enzyme_KSEA_diffexp_sub_FC/", 
                                             cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
        ksea_diffexp <- ksea_diffexp[ksea_diffexp$GENE %in% ksea_diffexp2$GENE,]
        rm(ksea_diffexp2)
        ## annotate regression results within tumours
        ksea_diffexp$reg.tumor <- "data_unavailable"
        tmp <- as.vector(ksea_diffexp$reg.tumor)
        tmp[!is.na(ksea_diffexp$FDR_pho_kin)] <- "other"
        tmp[!is.na(ksea_diffexp$FDR_pho_kin) & ksea_diffexp$enzyme_type == "kinase" & ksea_diffexp$P_pho_kin < reg_sig & ksea_diffexp$coef_pho_kin > 0] <- paste0("positive(p.value<",reg_sig, ")")
        tmp[!is.na(ksea_diffexp$FDR_pho_kin) & ksea_diffexp$enzyme_type == "kinase" & ksea_diffexp$FDR_pho_kin < reg_sig & ksea_diffexp$coef_pho_kin > 0] <- paste0("positive(FDR<",reg_sig, ")")
        tmp[!is.na(ksea_diffexp$FDR_pho_kin) & ksea_diffexp$enzyme_type == "phosphotase" & ksea_diffexp$P_pho_kin < reg_sig & ksea_diffexp$coef_pho_kin < 0] <- paste0("negative(p.value<",reg_sig, ")")
        tmp[!is.na(ksea_diffexp$FDR_pho_kin) & ksea_diffexp$enzyme_type == "phosphotase" & ksea_diffexp$FDR_pho_kin < reg_sig & ksea_diffexp$coef_pho_kin < 0] <- paste0("negative(FDR<",reg_sig, ")")
        ksea_diffexp$reg.tumor <- tmp
        
        ## annotate regression results within tumours+normals
        ksea_diffexp$reg.tn <- "data_unavailable"
        tmp <- as.vector(ksea_diffexp$reg.tn)
        tmp[!is.na(ksea_diffexp$FDR_pho_kin.tn)] <- "other"
        tmp[!is.na(ksea_diffexp$FDR_pho_kin.tn) & ksea_diffexp$enzyme_type == "kinase" & ksea_diffexp$P_pho_kin.tn < reg_sig & ksea_diffexp$coef_pho_kin.tn > 0] <- paste0("positive(p.value<",reg_sig, ")")
        tmp[!is.na(ksea_diffexp$FDR_pho_kin.tn) & ksea_diffexp$enzyme_type == "kinase" & ksea_diffexp$FDR_pho_kin.tn < reg_sig & ksea_diffexp$coef_pho_kin.tn > 0] <- paste0("positive(FDR<",reg_sig, ")")
        tmp[!is.na(ksea_diffexp$FDR_pho_kin.tn) & ksea_diffexp$enzyme_type == "phosphotase" & ksea_diffexp$P_pho_kin.tn < reg_sig & ksea_diffexp$coef_pho_kin.tn < 0] <- paste0("negative(p.value<",reg_sig, ")")
        tmp[!is.na(ksea_diffexp$FDR_pho_kin.tn) & ksea_diffexp$enzyme_type == "phosphotase" & ksea_diffexp$FDR_pho_kin.tn < reg_sig & ksea_diffexp$coef_pho_kin.tn < 0] <- paste0("negative(FDR<",reg_sig, ")")
        ksea_diffexp$reg.tn <- tmp
          
        ## annotate regression results within tumours+normals
        ksea_diffexp_direction <- data.frame(table(ksea_diffexp[, c("GENE", "substrate_FC_diffexp")]))
        
        ## get a list of enzymes with KSEA or differential expression data
        df <- unique(ksea_diffexp[!is.na(ksea_diffexp$enzyme_log2FC), c("GENE", "method", colnames(ksea_diffexp)[grepl(pattern = "enzyme", x = colnames(ksea_diffexp))])])
        df <- merge(df, ksea_diffexp_direction[, c("GENE", "substrate_FC_diffexp")], by = c("GENE"), all.x = T)

        ## annotate the combination of enzyme and substrate directions to look into
        df$enzyme_sub_consistent <- FALSE
        tmp <- as.vector(df$enzyme_sub_consistent)
        tmp[df$enzyme_type == "kinase" & df$enzyme_direction == "up" & grepl(pattern = "up", df$substrate_FC_diffexp)] <- TRUE
        tmp[df$enzyme_type == "kinase" & df$enzyme_direction == "down" & grepl(pattern = "down", df$substrate_FC_diffexp)] <- TRUE
        tmp[df$enzyme_type == "phosphotase" & df$method == "KSEA" & df$enzyme_log2FC > 0 & grepl(pattern = "up", df$substrate_FC_diffexp)] <- TRUE
        tmp[df$enzyme_type == "phosphotase" & df$method == "KSEA" & df$enzyme_log2FC < 0 & grepl(pattern = "down", df$substrate_FC_diffexp)] <- TRUE
        tmp[df$enzyme_type == "phosphotase" & df$method == "diffexp" & df$enzyme_log2FC < 0 & grepl(pattern = "up", df$substrate_FC_diffexp)] <- TRUE
        tmp[df$enzyme_type == "phosphotase" & df$method == "diffexp" & df$enzyme_log2FC > 0 & grepl(pattern = "down", df$substrate_FC_diffexp)] <- TRUE
        df$enzyme_sub_consistent <- tmp
        df <- df[df$enzyme_sub_consistent == TRUE,]
        
        ## annotate regression results within tumours+normals
        ksea_diffexp_reg <- data.frame(table(ksea_diffexp[!is.na(ksea_diffexp$enzyme_log2FC), c("GENE", "substrate_FC_diffexp", "reg.tumor")]))
        df <- merge(df, ksea_diffexp_reg, by = c("GENE", "substrate_FC_diffexp"), all.x = T)
        colnames(df)[1] <- "Kinase.Gene" 
        
        df$Freq2p <- df$Freq
        df$Freq2p[grepl(x = df$enzyme_direction, pattern = "down")] <- (-df$Freq2p[grepl(x = df$enzyme_direction, pattern = "down")])
        df$Freq2p <- as.numeric(df$Freq2p)
        df$substrate_direction <- ifelse(grepl(pattern = "up", x = df$substrate_FC_diffexp), "up", "down")
        x <- sapply(unique(df$Kinase.Gene), FUN = function(g, dat) {
          rows <- which(dat$Kinase.Gene == g)
          direction_g <- unique(dat$enzyme_direction[rows])
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
        df <- unique(df)
        
        for (enzyme_type in c("kinase")) {
          df_e <- df[df$Kinase.Gene %in% op_enzyme$GENE[op_enzyme$enzyme_type == enzyme_type],]
          df_e$Kinase.Gene <- factor(df_e$Kinase.Gene, levels = as.vector(df_e$Kinase.Gene)[order(as.vector(abs(df_e$enzyme_log2FC)))])
          df_e$enzyme_direction <- factor(df_e$enzyme_direction, levels = c("up", "down"))
          
          df_e_kin <- unique(df_e[, c("Kinase.Gene", "method", colnames(df_e)[grepl(pattern = "enzyme", x = colnames(df_e))])])
          df_e_kin$method <- factor(df_e_kin$method, levels = c("diffexp", "KSEA"))
          
          for (enzyme_direction in c("up", "down")) {
            subdir1 <- paste0(makeOutDir(), "facet_enzyme_direction/")
            dir.create(subdir1)
            subdir2 <- paste0(subdir1, "w.substrates/")
            dir.create(subdir2)
            df_e_d <- df_e[df_e$enzyme_direction == enzyme_direction,]
            df_e_d$Kinase.Gene <- factor(df_e_d$Kinase.Gene, levels = as.vector(df_e_d$Kinase.Gene)[order(as.vector(abs(df_e_d$enzyme_log2FC)))])
            
            # plot number of substrate phosphosites for each up/down-regulated kinases
            p <- ggplot()
            p <- p + geom_bar(data=df_e_d, aes(y = Freq2p, x = Kinase.Gene, fill = reg.tumor ), stat="identity")
            # p <- p + geom_text(data=df_e_d, aes(y = x, x = Kinase.Gene, label = Kinase.Gene, color = enzyme_direction), size = 2)
            p <- p + scale_color_manual(values = c("down" = "#1F78B4", "up" = "#E31A1C"))
            p <- p + scale_fill_manual(values = c("data_unavailable" = "#d9d9d9", "other" = "#969696", 
                                                  "positive(p.value<0.05)" = "#FB9A99", "positive(FDR<0.05)" = "#E31A1C",
                                                  "negative(p.value<0.05)" = "#A6CEE3", "negative(FDR<0.05)" = "#1F78B4"))
            p <- p + theme_bw() + theme_nogrid()
            p <- p + facet_grid(enzyme_direction~., scales = "free", space = "free_y")
            p <- p + coord_flip()
            p <- p + xlab("enzyme")+ylab("number of differentially phosphorylated phosphosites")
            p <- p + theme(axis.title=element_text(size=10))
            p <- p + theme(axis.text.x = element_text(colour="black", size=10),
                           axis.text.y = element_text(colour= ifelse(enzyme_direction == "up","#E31A1C", "#1F78B4"), size = 5))
            p
            fn = paste0(subdir2, "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "_", cancer, "_", enzyme_type, "_", enzyme_direction, "_KSEA_diffexp_substrates.pdf")
            ggsave(file=fn, height=6, width=5)
            
          }
        }
        
        for (enzyme_type in c("phosphotase")) {
          df_e <- df[df$Kinase.Gene %in% op_enzyme$GENE[op_enzyme$enzyme_type == enzyme_type],]
          df_e$enzyme_log2FC[df_e$method == "KSEA"] <- -(df_e$enzyme_log2FC[df_e$method == "KSEA"])
          df_e$Kinase.Gene <- factor(df_e$Kinase.Gene, levels = as.vector(df_e$Kinase.Gene)[order(as.vector(abs(df_e$enzyme_log2FC)))])
          df_e$enzyme_direction <- factor(df_e$enzyme_direction, levels = c("up", "down"))
          
          df_e_kin <- unique(df_e[, c("Kinase.Gene", "method", colnames(df_e)[grepl(pattern = "enzyme", x = colnames(df_e))])])
          df_e_kin$method <- factor(df_e_kin$method, levels = c("diffexp", "KSEA"))
          
          p <- ggplot()
          p <- p + geom_bar(data=df_e_kin, aes(y = enzyme_log2FC, x = Kinase.Gene, fill =  method_enzyme_direction), stat="identity")
          p <- p + geom_text(data=df_e_kin, aes(y = 0.6*enzyme_log2FC, x = Kinase.Gene, label = Kinase.Gene), size = 2, color = "black")
          p <- p + scale_fill_manual(values = c("diffexp-up" = "#FF7F00", "KSEA-up" = "#1F78B4", "diffexp-down" = "#CAB2D6", "KSEA-down"="#E31A1C"))
          p <- p + theme_bw() + theme_nogrid()
          p <- p + facet_grid(method_enzyme_direction~., scales = "free", space = "free_y")
          p <- p + coord_flip()
          p <- p + xlab(enzyme_type)+ylab("log2FC(tumor/normal)")
          p <- p + theme(axis.title=element_text(size=10))
          p <- p + theme(axis.text.x = element_text(colour="black", size=10),
                         axis.text.y = element_blank())
          p
          subdir <- paste0(makeOutDir(), "facet_enzyme_direction_method/")
          dir.create(subdir)
          fn = paste0(subdir, "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "_", cancer, "_", enzyme_type, "_KSEA_diffexp.pdf")
          
          h = 3
          ggsave(file=fn, height=h, width=4)
        }
      }
    }
  }
}
