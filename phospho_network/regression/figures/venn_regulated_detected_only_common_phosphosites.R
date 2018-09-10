# Yige Wu @ WashU 2018 Apr
# venn diagram showing whether regulated enzyme-substrate pairs only in breast are detected in colorectal and ovarian cancer

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(eulerr)


# set variables -----------------------------------------------------------
# reg_nonNA2test <- seq(from = 5, to = 30, by = 5)
reg_nonNA2test <- c(5, 25)
num_nonNA_common2test <- seq(from = 5, to = 30, by = 5)

# inputs ------------------------------------------------------------------


# plot venn diagram -------------------------------------------------------
regulated_ratio_across_thres <- NULL
for (num_nonNA_common in num_nonNA_common2test) {
  for (reg_nonNA in reg_nonNA2test) {
    for (enzyme_type in c("kinase")) {
      pho_sites_cans <- list()
      pho_sites_all <- NULL
      
      for (cancer in cancers_sort) {
        pho_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
        samples <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
        pho_data <- pho_data[rowSums(!is.na(pho_data[, samples])) >= num_nonNA_common,]
        pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
        pho_sites <- unique(pho_head[, c("SUBSTRATE", "SUB_MOD_RSD")])
        pho_sites$site <- paste0(pho_sites$SUBSTRATE, "_", pho_sites$SUB_MOD_RSD)
        pho_sites_cans[[cancer]] <- unique(as.vector(pho_sites$site))
        pho_sites_all <- unique(c(pho_sites_all, as.vector(pho_sites$site)))
      }
      dat <- data.frame(site = pho_sites_all)
      for (cancer in cancers_sort) {
        dat[, cancer] <- FALSE
        dat[dat$site %in% pho_sites_cans[[cancer]], cancer]  <- TRUE
      }
      pho_site_common <- dat[dat$BRCA & dat$OV & dat$CO,]
      
      sup_cans_tab_en <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/", 
                                              enzyme_type, "_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
      sup_cans_tab_en <- sup_cans_tab_en[sup_cans_tab_en$Size >= reg_nonNA,]
      name = c("pro_kin","pro_sub","pho_kin")
      
      ## adjust p-values to FDR
      for (cancer in cancers_sort) {
        for(self in c(TRUE,FALSE)) {
          for(coln in name) {#adjust pvalues for each variable
            row <- (sup_cans_tab_en$self==self) & (sup_cans_tab_en$Cancer==cancer)
            sup_cans_tab_en[row,paste("FDR_",coln,sep = "")] <-p.adjust(sup_cans_tab_en[row,paste("P_",coln,sep = "")],method = "fdr")
          }
        }
      }
      
      sup_cans_tab_en$GENE <- as.vector(sup_cans_tab_en$KINASE)
      sup_cans_tab_en$SUB_GENE <- as.vector(sup_cans_tab_en$SUBSTRATE)
      sup_cans_tab_en <- markSigKS(regression = sup_cans_tab_en, sig_thres = reg_sig[enzyme_type], enzyme_type = enzyme_type)
      sup_cans_tab_en$regulated <- (sup_cans_tab_en$fdr_sig & sup_cans_tab_en$coef_sig)
      
      ## filter by phosphosite
      sup_cans_tab_en$site <- paste0(sup_cans_tab_en$SUB_GENE, "_", sup_cans_tab_en$SUB_MOD_RSD)
      sup_cans_tab_en <- sup_cans_tab_en[sup_cans_tab_en$site %in% pho_site_common$site,]
      detected_cans <- NULL
      regulated_cans <- NULL
      pairs_all <- NULL
      detected_cans <- list()
      regulated_cans <- list()
      for (cancer in cancers_sort) {
        detected_cans[[cancer]] <- as.vector(unique(sup_cans_tab_en$pair[sup_cans_tab_en$Cancer == cancer]))
        regulated_cans[[cancer]] <- as.vector(unique(sup_cans_tab_en$pair[sup_cans_tab_en$regulated & sup_cans_tab_en$Cancer == cancer]))
        pairs_all <- unique(c(pairs_all, detected_cans[[cancer]]))
      }
      
      dat <- data.frame(pair = pairs_all)
      dat <- merge(dat, unique(sup_cans_tab_en[, c("pair", "SELF")]), all.x = T)
      for (cancer in cancers_sort) {
        dat[, paste0("detected_", cancer)] <- FALSE
        dat[dat$pair %in% detected_cans[[cancer]], paste0("detected_", cancer)]  <- TRUE
        
        dat[, paste0("regulated_", cancer)] <- FALSE
        dat[dat$pair %in% regulated_cans[[cancer]], paste0("regulated_", cancer)]  <- TRUE
      }
      dat <- data.frame(dat)
      regulated_stats <- data.frame(table(dat[,-1]))
      SELF <- NULL
      Cancer <- NULL
      count <- 0
      regulated_num <- NULL
      detected_num <- NULL
      regulated_ratio <- NULL
      for (self in c("cis", "trans")) {
        for (cancer2 in cancers_sort) {
          count <- count + 1
          SELF[count] <- self
          Cancer[count] <- cancer2
          regulated_num[count] <- sum(regulated_stats$Freq[regulated_stats$SELF == self & regulated_stats[, paste0("regulated_", cancer2)] == "TRUE"])
          detected_num[count] <- sum(regulated_stats$Freq[regulated_stats$SELF == self & regulated_stats[, paste0("detected_", cancer2)] == "TRUE"])
          regulated_ratio[count] <- regulated_num[count]/detected_num[count]
        }
      }
      regulated_ratio_per_thres <- data.frame(SELF, Cancer, regulated_num, detected_num, regulated_ratio)
      regulated_ratio_per_thres$reg_nonNA <- reg_nonNA
      regulated_ratio_per_thres$num_nonNA_common <- num_nonNA_common
      regulated_ratio_across_thres <- rbind(regulated_ratio_across_thres, regulated_ratio_per_thres)
      
      fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_substrate_pairs_regulated_commonSiteNonNA', num_nonNA_common, '_regresssionNonNA', reg_nonNA ,'_venn_cistrans.pdf',sep ="")
      grid.newpage()
      pdf(fn, height = 12, width = 20, useDingbats = FALSE)
      fit <- euler(combinations = dat[, c(1,2,4,6,8)], input = "disjoint", shape = 'circle', by = SELF)
      p <-plot(fit, quantities = list(fontsize = 30), fills = color_cancers2, legend = list(fontsize = 30))
      grid.draw(p)
      dev.off()
      
      # for (cancer in cancers_sort) {
      #   dat2p <- dat[, c("pair", paste0("detected_", cancers_sort), paste0("regulated_", cancer))]
      #   fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_substrate_pairs_detected_regualted_venn_', cancer, '.pdf',sep ="")
      #   grid.newpage()
      #   pdf(fn, height = 12, width = 12, useDingbats = FALSE)
      #   fit <- euler(combinations = dat2p, input = "disjoint", shape = 'circle')
      #   plot(fit, quantities = TRUE)
      #   p <-plot(fit, quantities = list(fontsize = 10), legend = list(fontsize = 30), labels = F, fills = c(color_cancers2, color_cancers2[cancer]))
      #   grid.draw(p)
      #   dev.off()
      # }
    }
  }
}


# changing the threshold of non-NA values ---------------------------------
## get list of detected and regulated
for (reg_nonNA in reg_nonNA2test) {
  for (self in c("cis")) {
    tab2p <- regulated_ratio_across_thres
    tab2p <- tab2p[tab2p$SELF == self & !is.na(tab2p$regulated_ratio) & tab2p$reg_nonNA == reg_nonNA,]
    p <- ggplot()
    p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA_common, y = detected_num/10, group = Cancer, color = Cancer), linetype = 2)
    p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA_common, y = regulated_ratio*100, group = Cancer, color = Cancer), linetype = 1)
    p <- p + scale_y_continuous(breaks = seq(from = 0, to = 150, by = 10), labels = seq(from = 0, to = 150, by = 10)*10,
                                sec.axis = sec_axis(~./100, name = "regulated [%]", breaks = seq(from = 0, to = 1, by = 0.1), labels = paste0(seq(from = 0, to = 1, by = 0.1)*100, "%")))
    p <- p + ylab("number of detected pairs")
    p <- p + scale_color_manual(values = color_cancers2)
    p <- p + theme_minimal()
    p
    fn <- paste0(makeOutDir(resultD = resultD), self, "_regulated_ratio_across_commonSiteNonNA_regresssionNonNA", reg_nonNA, "_wdetected.pdf")
    ggsave(filename = fn, width = 6, height = 5)
    
    p <- ggplot()
    p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA_common, y = regulated_num/10, group = Cancer, color = Cancer), linetype = 2)
    p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA_common, y = regulated_ratio*100, group = Cancer, color = Cancer), linetype = 1)
    p <- p + scale_y_continuous(breaks = seq(from = 0, to = 150, by = 10), labels = seq(from = 0, to = 150, by = 10)*10,
                                sec.axis = sec_axis(~./100, name = "regulated [%]", breaks = seq(from = 0, to = 1, by = 0.1), labels = paste0(seq(from = 0, to = 1, by = 0.1)*100, "%")))
    p <- p + ylab("number of regulated pairs")
    p <- p + scale_color_manual(values = color_cancers2)
    p <- p + theme_minimal()
    p
    fn <- paste0(makeOutDir(resultD = resultD), self, "_regulated_ratio_across_thres_commonSiteNonNA_regresssionNonNA", reg_nonNA, "_wregulated.pdf")
    ggsave(filename = fn, width = 6, height = 5)
  }
  
  for (self in c("trans")) {
    tab2p <- regulated_ratio_across_thres
    tab2p <- tab2p[tab2p$SELF == self & !is.na(tab2p$regulated_ratio) & tab2p$reg_nonNA == reg_nonNA,]
    
    p <- ggplot()
    # p <- p + geom_point(data = tab2p, mapping = aes(x = reg_nonNA, y = y, group = Cancer, color = Cancer))
    p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA_common, y = detected_num/4000, group = Cancer, color = Cancer), linetype = 2)
    p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA_common, y = regulated_ratio*1000, group = Cancer, color = Cancer), linetype = 1)
    p <- p + scale_y_continuous(breaks = seq(from = 0, to = 60, by = 10), labels = seq(from = 0, to = 60, by = 10)*4000,
                                sec.axis = sec_axis(~./1000, name = "regulated [%]", breaks = seq(from = 0, to = 0.1, by = 0.01), labels = paste0(seq(from = 0, to = 0.1, by = 0.01)*100, "%")))
    p <- p + ylab("number of detected pairs")
    p <- p + scale_color_manual(values = color_cancers2)
    p <- p + theme_minimal()
    p
    fn <- paste0(makeOutDir(resultD = resultD), self, "_regulated_ratio_across_thres_commonSiteNonNA_regresssionNonNA", reg_nonNA, "_wdetected.pdf")
    ggsave(filename = fn, width = 6, height = 5)
    
    p <- ggplot()
    p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA_common, y = regulated_num/50, group = Cancer, color = Cancer), linetype = 2)
    p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA_common, y = regulated_ratio*500, group = Cancer, color = Cancer), linetype = 1)
    p <- p + scale_y_continuous(breaks = seq(from = 0, to = 60, by = 10), labels = seq(from = 0, to = 60, by = 10)*50,
                                sec.axis = sec_axis(~./500, name = "regulated [%]", breaks = seq(from = 0, to = 0.1, by = 0.01), labels = paste0(seq(from = 0, to = 0.1, by = 0.01)*100, "%")))
    p <- p + ylab("number of regulated pairs")
    p <- p + scale_color_manual(values = color_cancers2)
    p <- p + theme_minimal()
    p
    fn <- paste0(makeOutDir(resultD = resultD), self, "_regulated_ratio_across_thres_commonSiteNonNA_regresssionNonNA", reg_nonNA, "_wregulated.pdf")
    ggsave(filename = fn, width = 6, height = 5)
  }
}
