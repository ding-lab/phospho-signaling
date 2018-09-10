# Yige Wu @ WashU 2018 Jul
# summarize the prevalence of up-regulated and down-regulated ks pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(eulerr)

# variables ---------------------------------------------------------------
cancers_sort <- c("BRCA", "OV", "CO")
reg_sig <- 0.05
diff_sig <- 0.2
diff_log2fc <- 1

# inputs -------------------------------------------------------------------
sup_cans_tab <- NULL
for (cancer in cancers_sort) {
  sup_tab <- fread(paste0(ppnD, "kinase_activity/tables/fisher_es_pairs/", 
                          cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  sup_tab$Cancer <- cancer
  sup_cans_tab <- rbind(sup_cans_tab, sup_tab)
}

# summarize kinase-substrate pairs ----------------------------------------
enzyme_type <- "kinase"
sup_cans_tab_pk <- sup_cans_tab[sup_cans_tab$enzyme_type == enzyme_type,]
sup_cans_tab_pk <- data.frame(sup_cans_tab_pk)
sup_cans_tab_pk <- markSigSiteCan(sup_cans_tab_pk, sig_thres = reg_sig, enzyme_type = enzyme_type)
enzyme_type <- "phosphotase"
sup_cans_tab_pp <- sup_cans_tab[sup_cans_tab$enzyme_type == enzyme_type,]
sup_cans_tab_pp <- data.frame(sup_cans_tab_pp)
sup_cans_tab_pp <- markSigSiteCan(sup_cans_tab_pp, sig_thres = reg_sig, enzyme_type = enzyme_type)
sup_cans_tab <- rbind(sup_cans_tab_pp, sup_cans_tab_pk)

## annotate kinase substrate regulation
sup_cans_tab$regulated <- (sup_cans_tab$coef_sig & sup_cans_tab$fdr_sig)
enzyme_type <- "kinase"
sup_cans_tab_pk <- sup_cans_tab[sup_cans_tab$enzyme_type == enzyme_type,]
enzyme_type <- "phosphotase"
sup_cans_tab_pp <- sup_cans_tab[sup_cans_tab$enzyme_type == enzyme_type,]

## annotate up and down regulation
sup_cans_tab_pk$diffexp_type <- "other"
diffexp_sig <- ((!is.na(sup_cans_tab_pk$KSEA_pvalue) & sup_cans_tab_pk$KSEA_pvalue < diff_sig) | (!is.na(sup_cans_tab_pk$diffexp_log2FC & abs(sup_cans_tab_pk$diffexp_log2FC) > diff_log2fc)))
sup_cans_tab_pk$diffexp_type[diffexp_sig & sup_cans_tab_pk$enzyme_direction == "up" & sup_cans_tab_pk$substrate_direction == "up"] <- "up"
sup_cans_tab_pk$diffexp_type[diffexp_sig & sup_cans_tab_pk$enzyme_direction == "down" & sup_cans_tab_pk$substrate_direction == "down"] <- "down"
sup_cans_tab_pp$diffexp_type <- "other"
diffexp_sig <- ((!is.na(sup_cans_tab_pp$KSEA_pvalue) & sup_cans_tab_pp$KSEA_pvalue < diff_sig) | (!is.na(sup_cans_tab_pp$diffexp_log2FC & abs(sup_cans_tab_pp$diffexp_log2FC) > diff_log2fc)))
sup_cans_tab_pp$diffexp_type[diffexp_sig & sup_cans_tab_pp$enzyme_direction == "up" & sup_cans_tab_pp$substrate_direction == "up"] <- "down"
sup_cans_tab_pp$diffexp_type[diffexp_sig & sup_cans_tab_pp$enzyme_direction == "down" & sup_cans_tab_pp$substrate_direction == "down"] <- "up"

## get table
sum_tab1 <- data.frame(table(unique(sup_cans_tab_pk[, c("pair", "regulated", "diffexp_type", "Cancer")])[, c("regulated", "diffexp_type", "Cancer")]))
sum_tab2 <- data.frame(kinase_substrate_pair = c("regulated", "up", "down"))
sum_tab3 <- data.frame(table(unique(sup_cans_tab_pk[, c("pair", "regulated", "Cancer")])[, c("regulated", "Cancer")]))

for (cancer in cancers_sort) {
  sum_tab2[, cancer] = c(sum(sum_tab1$Freq[sum_tab1$regulated == "TRUE" & sum_tab1$Cancer == cancer]),
                         sum(sum_tab1$Freq[sum_tab1$diffexp_type == "up"& sum_tab1$Cancer == cancer]),
                         sum(sum_tab1$Freq[sum_tab1$diffexp_type == "down" & sum_tab1$Cancer == cancer]))
}
write.table(x = sum_tab2, file = paste0(makeOutDir(resultD = resultD), enzyme_type, "_substrate_pairs_table.txt"), 
            quote = F, row.names = F, sep = '\t')


# plot venn diagram per cancer, regulated & up & down-------------------------------------------------------
enzyme_type <- "kinase"
combs <- NULL
for (cancer in cancers_sort) {
  fn = paste(makeOutDir(resultD = resultD), cancer, '_', enzyme_type, '_substrate_pair_class_venn','.pdf',sep ="")
  grid.newpage()
  pdf(fn, height = 6, width = 6, useDingbats = FALSE)
  fit <- euler(combinations = c('regulated' = sum_tab1$Freq[sum_tab1$regulated == "TRUE" & sum_tab1$diffexp_type == "other" & sum_tab1$Cancer == cancer], 
                                'up' = sum_tab1$Freq[sum_tab1$regulated == "FALSE" & sum_tab1$diffexp_type == "up" & sum_tab1$Cancer == cancer], 
                                'down' = sum_tab1$Freq[sum_tab1$regulated == "FALSE" & sum_tab1$diffexp_type == "down" & sum_tab1$Cancer == cancer],
                                "regulated&up" = sum_tab1$Freq[sum_tab1$regulated == "TRUE" & sum_tab1$diffexp_type == "up" & sum_tab1$Cancer == cancer],
                                "regulated&down" = sum_tab1$Freq[sum_tab1$regulated == "TRUE" & sum_tab1$diffexp_type == "down" & sum_tab1$Cancer == cancer]), 
               input = "disjoint", shape = 'circle')
  p <-plot(fit, quantities = TRUE, fills = c('#CAB2D6', '#FB9A99', '#A6CEE3'))
  grid.draw(p)
  dev.off()
  
  ## try to get 3 cancers to plotted tgt
  tmp <- c(sum_tab1$Freq[sum_tab1$regulated == "TRUE" & sum_tab1$diffexp_type == "other" & sum_tab1$Cancer == cancer], 
           sum_tab1$Freq[sum_tab1$regulated == "FALSE" & sum_tab1$diffexp_type == "up" & sum_tab1$Cancer == cancer], 
           sum_tab1$Freq[sum_tab1$regulated == "FALSE" & sum_tab1$diffexp_type == "down" & sum_tab1$Cancer == cancer],
           sum_tab1$Freq[sum_tab1$regulated == "TRUE" & sum_tab1$diffexp_type == "up" & sum_tab1$Cancer == cancer],
           sum_tab1$Freq[sum_tab1$regulated == "TRUE" & sum_tab1$diffexp_type == "down" & sum_tab1$Cancer == cancer])
  names(tmp) <- c(paste0(cancer,'-regulated'), 
                  paste0(cancer,'-up'),
                  paste0(cancer,'-down'),
                  paste0(cancer,'-regulated&', cancer, '-up'),
                  paste0(cancer,'-regulated&', cancer, '-down'))
  combs <- c(combs, tmp)
}
fn = paste(makeOutDir(resultD = resultD), 'cancers', '_', enzyme_type, '_substrate_pair_class_venn','.pdf',sep ="")
grid.newpage()
pdf(fn, height = 12, width = 12, useDingbats = FALSE)
fit <- euler(combinations = combs, 
             input = "disjoint", shape = 'circle')

p <-plot(fit, quantities = TRUE, fills = c('#CAB2D6', '#FB9A99', '#A6CEE3', '#CAB2D6', '#FB9A99', '#A6CEE3', '#CAB2D6', '#FB9A99', '#A6CEE3'))
grid.draw(p)
dev.off()


# plot venn across cancers, only regulated --------------------------------
fn = paste(makeOutDir(resultD = resultD), 'cancers', '_', enzyme_type, '_substrate_pair_regulated_venn','.pdf',sep ="")
grid.newpage()
pdf(fn, height = 12, width = 12, useDingbats = FALSE)
all_pairs <- unique(sup_cans_tab_pk$pair[sup_cans_tab_pk$regulated & !is.na(sup_cans_tab_pk$regulated)])
dat <- data.frame(BRCA_regulated = (all_pairs %in% sup_cans_tab_pk$pair[sup_cans_tab_pk$sig_BRCA]),
                  OV_regulated = (all_pairs %in% sup_cans_tab_pk$pair[sup_cans_tab_pk$sig_OV]),
                  CO_regulated = (all_pairs %in% sup_cans_tab_pk$pair[sup_cans_tab_pk$sig_CO]))

dat_euler <- euler(dat)
p <- plot(dat_euler, quantities = TRUE, fills = c(colors['BRCA'], colors['OV'], colors['COAD']))
grid.draw(p)
dev.off()

# plot venn across cancers, up or down regulated --------------------------------
for (diffexp_type in c("up")) {
  fn = paste(makeOutDir(resultD = resultD), 'cancers', '_', enzyme_type, '_substrate_pair_', diffexp_type,'_venn','.pdf',sep ="")
  grid.newpage()
  pdf(fn, height = 12, width = 12, useDingbats = FALSE)
  all_pairs <- unique(sup_cans_tab_pk$pair)
  dat <- data.frame(BRCA_up = (all_pairs %in% sup_cans_tab_pk$pair[sup_cans_tab_pk$diffexp_type == diffexp_type & sup_cans_tab_pk$Cancer == "BRCA"]),
                    OV_up = (all_pairs %in% sup_cans_tab_pk$pair[sup_cans_tab_pk$diffexp_type == diffexp_type & sup_cans_tab_pk$Cancer == "OV"]),
                    CO_up = (all_pairs %in% sup_cans_tab_pk$pair[sup_cans_tab_pk$diffexp_type == diffexp_type & sup_cans_tab_pk$Cancer == "CO"]))
  
  dat_euler <- euler(dat)
  p <- plot(dat_euler, quantities = TRUE, fills = c(colors['BRCA'], colors['OV'], colors['COAD']))
  grid.draw(p)
  dev.off()
}

for (diffexp_type in c("down")) {
  fn = paste(makeOutDir(resultD = resultD), 'cancers', '_', enzyme_type, '_substrate_pair_', diffexp_type,'_venn','.pdf',sep ="")
  grid.newpage()
  pdf(fn, height = 12, width = 12, useDingbats = FALSE)
  all_pairs <- unique(sup_cans_tab_pk$pair)
  dat <- data.frame(BRCA_down = (all_pairs %in% sup_cans_tab_pk$pair[sup_cans_tab_pk$diffexp_type == diffexp_type & sup_cans_tab_pk$Cancer == "BRCA"]),
                    OV_down = (all_pairs %in% sup_cans_tab_pk$pair[sup_cans_tab_pk$diffexp_type == diffexp_type & sup_cans_tab_pk$Cancer == "OV"]),
                    CO_down = (all_pairs %in% sup_cans_tab_pk$pair[sup_cans_tab_pk$diffexp_type == diffexp_type & sup_cans_tab_pk$Cancer == "CO"]))
  
  dat_euler <- euler(dat)
  p <- plot(dat_euler, quantities = TRUE, fills = c(colors['BRCA'], colors['OV'], colors['COAD']))
  grid.draw(p)
  dev.off()
}
