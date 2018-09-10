# Yige Wu @ WashU 2018 Apr
# venn diagram showing whether regulated enzyme-substrate pairs only in breast are detected in colorectal and ovarian cancer

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(eulerr)


# set variables -----------------------------------------------------------

# inputs ------------------------------------------------------------------
## input enzyme-substrate pairs examined
regulated_ratio_across_thres <- NULL
for (sample_thres in seq(from = 5, to = 120, by = 5)) {
  for (enzyme_type in c("kinase")) {
    sup_cans_tab_en <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/", 
                                            enzyme_type, "_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
    sup_cans_tab_en <- sup_cans_tab_en[sup_cans_tab_en$Size >= sample_thres,]
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
    regulated_ratio_per_thres$sample_thres <- sample_thres
    regulated_ratio_across_thres <- rbind(regulated_ratio_across_thres, regulated_ratio_per_thres)
    
    fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_substrate_pairs_regulated_samples', sample_thres ,'_venn_cistrans.pdf',sep ="")
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


# changing the threshold of non-NA values ---------------------------------
## get list of detected and regulated
for (self in c("cis")) {
  tab2p <- regulated_ratio_across_thres[regulated_ratio_across_thres$SELF == self & !is.na(regulated_ratio_across_thres$regulated_ratio),]
  p <- ggplot()
  p <- p + geom_line(data = tab2p, mapping = aes(x = sample_thres, y = detected_num/10, group = Cancer, color = Cancer), linetype = 2)
  p <- p + geom_line(data = tab2p, mapping = aes(x = sample_thres, y = regulated_ratio*100, group = Cancer, color = Cancer), linetype = 1)
  p <- p + scale_y_continuous(breaks = seq(from = 0, to = 150, by = 10), labels = seq(from = 0, to = 150, by = 10)*10,
                              sec.axis = sec_axis(~./100, name = "regulated [%]", breaks = seq(from = 0, to = 1, by = 0.1), labels = paste0(seq(from = 0, to = 1, by = 0.1)*100, "%")))
  p <- p + ylab("number of detected pairs")
  p <- p + scale_color_manual(values = color_cancers2)
  p <- p + theme_minimal()
  p
  fn <- paste0(makeOutDir(resultD = resultD), self, "_regulated_ratio_across_thres_wdetected.pdf")
  ggsave(filename = fn, width = 6, height = 5)
  
  p <- ggplot()
  p <- p + geom_line(data = tab2p, mapping = aes(x = sample_thres, y = regulated_num/10, group = Cancer, color = Cancer), linetype = 1)
  p <- p + geom_line(data = tab2p, mapping = aes(x = sample_thres, y = regulated_ratio*100, group = Cancer, color = Cancer), linetype = 1)
  p <- p + scale_y_continuous(breaks = seq(from = 0, to = 150, by = 10), labels = seq(from = 0, to = 150, by = 10)*10,
                              sec.axis = sec_axis(~./100, name = "regulated [%]", breaks = seq(from = 0, to = 1, by = 0.1), labels = paste0(seq(from = 0, to = 1, by = 0.1)*100, "%")))
  p <- p + ylab("number of regulated pairs")
  p <- p + scale_color_manual(values = color_cancers2)
  p <- p + theme_minimal()
  p
  fn <- paste0(makeOutDir(resultD = resultD), self, "_regulated_ratio_across_thres_wregulated.pdf")
  ggsave(filename = fn, width = 6, height = 5)
}

for (self in c("trans")) {
  tab2p <- regulated_ratio_across_thres[regulated_ratio_across_thres$SELF == self & !is.na(regulated_ratio_across_thres$regulated_ratio),]
  
  p <- ggplot()
  # p <- p + geom_point(data = tab2p, mapping = aes(x = sample_thres, y = y, group = Cancer, color = Cancer))
  p <- p + geom_line(data = tab2p, mapping = aes(x = sample_thres, y = detected_num/4000, group = Cancer, color = Cancer), linetype = 2)
  p <- p + geom_line(data = tab2p, mapping = aes(x = sample_thres, y = regulated_ratio*1000, group = Cancer, color = Cancer), linetype = 1)
  p <- p + scale_y_continuous(breaks = seq(from = 0, to = 60, by = 10), labels = seq(from = 0, to = 60, by = 10)*4000,
                              sec.axis = sec_axis(~./1000, name = "regulated [%]", breaks = seq(from = 0, to = 0.1, by = 0.01), labels = paste0(seq(from = 0, to = 0.1, by = 0.01)*100, "%")))
  p <- p + ylab("number of detected pairs")
  p <- p + scale_color_manual(values = color_cancers2)
  p <- p + theme_minimal()
  p
  fn <- paste0(makeOutDir(resultD = resultD), self, "_regulated_ratio_across_thres_wdetected_num.pdf")
  ggsave(filename = fn, width = 6, height = 5)
  
  p <- ggplot()
  p <- p + geom_line(data = tab2p, mapping = aes(x = sample_thres, y = regulated_num/50, group = Cancer, color = Cancer), linetype = 1)
  p <- p + geom_line(data = tab2p, mapping = aes(x = sample_thres, y = regulated_ratio*1000, group = Cancer, color = Cancer), linetype = 1)
  p <- p + scale_y_continuous(breaks = seq(from = 0, to = 60, by = 10), labels = seq(from = 0, to = 60, by = 10)*50,
                              sec.axis = sec_axis(~./1000, name = "regulated [%]", breaks = seq(from = 0, to = 0.1, by = 0.01), labels = paste0(seq(from = 0, to = 0.1, by = 0.01)*100, "%")))
  p <- p + ylab("number of regulated pairs")
  p <- p + scale_color_manual(values = color_cancers2)
  p <- p + theme_minimal()
  p
  fn <- paste0(makeOutDir(resultD = resultD), self, "_regulated_ratio_across_thres_wregulated_num.pdf")
  ggsave(filename = fn, width = 6, height = 5)
}