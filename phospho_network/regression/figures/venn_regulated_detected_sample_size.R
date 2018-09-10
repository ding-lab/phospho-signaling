# Yige Wu @ WashU 2018 Apr
# venn diagram showing whether regulated enzyme-substrate pairs only in breast are detected in colorectal and ovarian cancer

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(eulerr)
library(ggrepel)
# set variables -----------------------------------------------------------
num_nonNA2test <- seq(from = 5, to = 120, by = 5)
# num_nonNA2test <- seq(from = 5, to = 120, by = 5)

# inputs ------------------------------------------------------------------
## input enzyme-substrate pairs examined
regulated_ratio_across_thres <- NULL
for (time in 1:1) {
  for (num_nonNA in num_nonNA2test) {
    for (enzyme_type in c("kinase")) {
      downsize_tab <- NULL
      for (cancer in c("BRCA", "CO")) {
        can_tab <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level_BRCA_downsize/", 
                                                time, "/", enzyme_type, "_substrate_regression_cptac2p_", cancer, "_tumor.txt"), data.table = F)
        downsize_tab <- rbind(can_tab, downsize_tab)
      }
      original_tab <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/", 
                                              enzyme_type, "_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
      downsize_tab <- rbind(downsize_tab, original_tab[original_tab$Cancer == "OV",])
      downsize_tab <- downsize_tab[downsize_tab$Size >= num_nonNA,]
      original_tab <- original_tab[original_tab$Size >= num_nonNA,]
      name = c("pro_kin","pro_sub","pho_kin")
      
      ## adjust p-values to FDR
      for (cancer in cancers_sort) {
        for(self in c(TRUE,FALSE)) {
          for(coln in name) {#adjust pvalues for each variable
            row <- (downsize_tab$self==self) & (downsize_tab$Cancer==cancer)
            downsize_tab[row,paste("FDR_",coln,sep = "")] <-p.adjust(downsize_tab[row,paste("P_",coln,sep = "")],method = "fdr")
          }
        }
      }
      
      downsize_tab$GENE <- as.vector(downsize_tab$KINASE)
      downsize_tab$SUB_GENE <- as.vector(downsize_tab$SUBSTRATE)
      downsize_tab <- markSigKS(regression = downsize_tab, sig_thres = reg_sig[enzyme_type], enzyme_type = enzyme_type)
      downsize_tab$regulated <- (downsize_tab$fdr_sig & downsize_tab$coef_sig)
      
      original_tab$GENE <- as.vector(original_tab$KINASE)
      original_tab$SUB_GENE <- as.vector(original_tab$SUBSTRATE)
      original_tab <- markSigKS(regression = original_tab, sig_thres = reg_sig[enzyme_type], enzyme_type = enzyme_type)
      original_tab$regulated <- (original_tab$fdr_sig & original_tab$coef_sig)
      
      ## get the pairs disappear after down-sizing
      original_reg <- original_tab[original_tab$regulated,]
      downsize_lost <- original_reg
      downsize_lost <- downsize_lost[!(paste0(downsize_lost$pair, downsize_lost$Cancer) %in% paste0(downsize_tab$pair[downsize_tab$regulated], downsize_tab$Cancer[downsize_tab$regulated])),]
      write.table(x = downsize_lost, file = paste0(makeOutDir(resultD = resultD),
                                                   enzyme_type, '_substrate_pairs_regulated_samples', num_nonNA ,'_downsizelost', time, ".txt"),
                  row.names = F, col.names = T, quote = F, sep = "\t")
      for (self in "cis") {
        tab_foreground <- downsize_lost[downsize_lost$SELF == self,]
        tab_foreground$log10_FDR <- -log10(tab_foreground$FDR_pro_kin)
        tab_foreground$coef <- as.vector(tab_foreground$coef_pro_kin)
        
        tab_backgroud <- original_reg[original_reg$SELF == self & !(original_reg$Cancer == "OV"),]
        tab_backgroud$log10_FDR <- -log10(tab_backgroud$FDR_pro_kin)
        tab_backgroud$coef <- as.vector(tab_backgroud$coef_pro_kin)
        tab_backgroud <- tab_backgroud[tab_backgroud$coef < max(tab_foreground$coef),]

        
        p <- ggplot()
        p <- p + geom_point(data = tab_backgroud, mapping = aes(x = coef, y = log10_FDR, size = Size/100), 
                            alpha = 0.2, color = "black", shape = 21, fill = "grey")
        p <- p + geom_point(data = tab_foreground, mapping = aes(x = coef, y = log10_FDR, size = Size/100, fill = Cancer), 
                            alpha = 0.5, color = "black", shape = 21)
        p <- p + scale_fill_manual(values = color_cancers2)
        p <- p + scale_color_manual(values = color_cancers2)
        p <- p + geom_text_repel(data = tab_foreground[tab_foreground$coef > quantile(x = tab_foreground$coef, probs = 0.99) | tab_foreground$log10_FDR > quantile(x = tab_foreground$log10_FDR, 0.99),], 
                                 mapping = aes(x = coef, y = log10_FDR, label = pair), size = 3, alpha = 1, force = 1)
        p <- p + facet_grid(.~Cancer, space = "fixed", scales = "fixed")
        p <- p + theme_nogrid()
        p <- p + ylim(c(-5, 50))
        p
        ggsave(filename = paste0(makeOutDir(resultD = resultD), self, "_",
                                 enzyme_type, '_substrate_pairs_regulated_samples', num_nonNA ,'_downsizelost', time, ".pdf"),
               width = 8, height = 4)
      }
      
      for (self in "trans") {
        tab_foreground <- downsize_lost[downsize_lost$SELF == self,]
        tab_foreground$log10_FDR <- -log10(tab_foreground$FDR_pho_kin)
        tab_foreground$coef <- as.vector(tab_foreground$coef_pho_kin)
        
        tab_backgroud <- original_reg[original_reg$SELF == self & !(original_reg$Cancer == "OV"),]
        tab_backgroud$log10_FDR <- -log10(tab_backgroud$FDR_pho_kin)
        tab_backgroud$coef <- as.vector(tab_backgroud$coef_pho_kin)
        tab_backgroud <- tab_backgroud[tab_backgroud$coef < max(tab_foreground$coef),]
        
        
        p <- ggplot()
        p <- p + geom_point(data = tab_backgroud, mapping = aes(x = coef, y = log10_FDR, size = Size/100), 
                            alpha = 0.2, color = "black", shape = 21, fill = "grey")
        p <- p + geom_point(data = tab_foreground, mapping = aes(x = coef, y = log10_FDR, size = Size/100, fill = Cancer), 
                            alpha = 0.5, color = "black", shape = 21)
        p <- p + scale_fill_manual(values = color_cancers2)
        p <- p + scale_color_manual(values = color_cancers2)
        p <- p + geom_text_repel(data = tab_foreground[tab_foreground$coef > quantile(x = tab_foreground$coef, probs = 0.995),],
                                 mapping = aes(x = coef, y = log10_FDR, label = pair), size = 3, alpha = 1, force = 1)
        p <- p + facet_grid(.~Cancer, space = "fixed", scales = "fixed")
        p <- p + theme_nogrid()
        p <- p + ylim(c(-1, 20))
        p
        ggsave(filename = paste0(makeOutDir(resultD = resultD), self, "_",
                                 enzyme_type, '_substrate_pairs_regulated_samples', num_nonNA ,'_downsizelost', time, ".pdf"),
               width = 8, height = 4)
      }
      
      
      detected_cans <- NULL
      regulated_cans <- NULL
      pairs_all <- NULL
      detected_cans <- list()
      regulated_cans <- list()
      for (cancer in cancers_sort) {
        detected_cans[[cancer]] <- as.vector(unique(downsize_tab$pair[downsize_tab$Cancer == cancer]))
        regulated_cans[[cancer]] <- as.vector(unique(downsize_tab$pair[downsize_tab$regulated & downsize_tab$Cancer == cancer]))
        pairs_all <- unique(c(pairs_all, detected_cans[[cancer]]))
      }
      
      dat <- data.frame(pair = pairs_all)
      dat <- merge(dat, unique(downsize_tab[, c("pair", "SELF")]), all.x = T)
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
      regulated_ratio_per_thres$num_nonNA <- num_nonNA
      regulated_ratio_across_thres <- rbind(regulated_ratio_across_thres, regulated_ratio_per_thres)
      
      fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_substrate_pairs_regulated_samples', num_nonNA ,'_venn_cistrans_downsize', time, '.pdf',sep ="")
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

## get list of detected and regulated
for (self in c("cis")) {
  tab2p <- regulated_ratio_across_thres[regulated_ratio_across_thres$SELF == self & !is.na(regulated_ratio_across_thres$regulated_ratio),]
  p <- ggplot()
  p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA, y = detected_num/10, group = Cancer, color = Cancer), linetype = 2)
  p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA, y = regulated_ratio*100, group = Cancer, color = Cancer), linetype = 1)
  p <- p + scale_y_continuous(breaks = seq(from = 0, to = 150, by = 10), labels = seq(from = 0, to = 150, by = 10)*10,
                              sec.axis = sec_axis(~./100, name = "regulated [%]", breaks = seq(from = 0, to = 1, by = 0.1), labels = paste0(seq(from = 0, to = 1, by = 0.1)*100, "%")))
  p <- p + ylab("number of detected pairs")
  p <- p + scale_color_manual(values = color_cancers2)
  p <- p + theme_minimal()
  p
  fn <- paste0(makeOutDir(resultD = resultD), self, "_regulated_ratio_across_thres_wdetected.pdf")
  ggsave(filename = fn, width = 6, height = 5)
  
  p <- ggplot()
  p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA, y = regulated_num/10, group = Cancer, color = Cancer), linetype = 1)
  p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA, y = regulated_ratio*100, group = Cancer, color = Cancer), linetype = 1)
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
  # p <- p + geom_point(data = tab2p, mapping = aes(x = num_nonNA, y = y, group = Cancer, color = Cancer))
  p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA, y = detected_num/4000, group = Cancer, color = Cancer), linetype = 2)
  p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA, y = regulated_ratio*1000, group = Cancer, color = Cancer), linetype = 1)
  p <- p + scale_y_continuous(breaks = seq(from = 0, to = 60, by = 10), labels = seq(from = 0, to = 60, by = 10)*4000,
                              sec.axis = sec_axis(~./1000, name = "regulated [%]", breaks = seq(from = 0, to = 0.1, by = 0.01), labels = paste0(seq(from = 0, to = 0.1, by = 0.01)*100, "%")))
  p <- p + ylab("number of detected pairs")
  p <- p + scale_color_manual(values = color_cancers2)
  p <- p + theme_minimal()
  p
  fn <- paste0(makeOutDir(resultD = resultD), self, "_regulated_ratio_across_thres_wdetected_num.pdf")
  ggsave(filename = fn, width = 6, height = 5)
  
  p <- ggplot()
  p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA, y = regulated_num/50, group = Cancer, color = Cancer), linetype = 1)
  p <- p + geom_line(data = tab2p, mapping = aes(x = num_nonNA, y = regulated_ratio*1000, group = Cancer, color = Cancer), linetype = 1)
  p <- p + scale_y_continuous(breaks = seq(from = 0, to = 60, by = 10), labels = seq(from = 0, to = 60, by = 10)*50,
                              sec.axis = sec_axis(~./1000, name = "regulated [%]", breaks = seq(from = 0, to = 0.1, by = 0.01), labels = paste0(seq(from = 0, to = 0.1, by = 0.01)*100, "%")))
  p <- p + ylab("number of regulated pairs")
  p <- p + scale_color_manual(values = color_cancers2)
  p <- p + theme_minimal()
  p
  fn <- paste0(makeOutDir(resultD = resultD), self, "_regulated_ratio_across_thres_wregulated_num.pdf")
  ggsave(filename = fn, width = 6, height = 5)
}