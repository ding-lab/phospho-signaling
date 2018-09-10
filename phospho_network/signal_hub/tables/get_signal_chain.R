cancer <- "CO"
ksea <- fread(input = paste0("/Users/yigewu/Box\ Sync/", "cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA/OmniPath_NetworKIN_", 
                             "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m.cutoff, "/", cancer, "/KSEA Kinase Scores.csv"),
              data.table = F)              
enzyme <- "CDK1"

up_tree <- vector("list")
kin_tmp <-  enzyme
up_tree[[kin_tmp]] <- vector("list")
subs_tmp <- unique(ptms_site_pairs_sup$SUB_GENE[ptms_site_pairs_sup$GENE == kin_tmp & ptms_site_pairs_sup$SUB_GENE != kin_tmp])
kins_tmp <- rep(c(kin_tmp), c(length(unique(ptms_site_pairs_sup$SUB_GENE[ptms_site_pairs_sup$GENE == kin_tmp]))))
subs_current_tmp <- subs_tmp
kins_current_tmp <- kins_tmp
count <- 0
while (length(subs_tmp) > 0 & count < 1000) {
  # print(subs_tmp)
  count <- count + 1
  subs_tmp <- subs_tmp[subs_tmp!=kins_tmp]
  kins_tmp <- kins_tmp[kins_tmp!=subs_tmp]
  while (kins_tmp[1] %in% names(up_tree) && subs_tmp[1] %in% names(up_tree[[kins_tmp[1]]])) {
    subs_tmp <- subs_tmp[-1]
    kins_tmp <- kins_tmp[-1]
  }
  sub_tmp <- subs_tmp[1]
  kin_tmp <- kins_tmp[1]
  ksea_sub <- ksea[ksea$Kinase.Gene == sub_tmp & ksea$p.value < 0.2,]
  pho_diff_sub <- Pho.diffexp3can[Pho.diffexp3can$Cancer == cancer & Pho.diffexp3can$Gene == sub_tmp & Pho.diffexp3can$`Fold Change` > 1,]
  if ((nrow(ksea_sub) > 0 && ksea_sub$z.score > 0) || nrow(pho_diff_sub) > 0) {
  # if (nrow(pho_diff_sub) > 0) {
    cat(paste0(kin_tmp, " -> ", sub_tmp, "\n"))
    up_tree[[kin_tmp]][[sub_tmp]][["ksea"]] <- ksea_sub
    up_tree[[kin_tmp]][[sub_tmp]][["diffexp"]] <- pho_diff_sub
    subs_current_tmp <- unique(ptms_site_pairs_sup$SUB_GENE[ptms_site_pairs_sup$GENE == sub_tmp & ptms_site_pairs_sup$SUB_GENE != sub_tmp])
    kins_current_tmp <- rep(c(sub_tmp), c(length(subs_current_tmp)))
    subs_tmp <- c(subs_current_tmp, subs_tmp)
    kins_tmp <- c(kins_current_tmp, kins_tmp)
    if (length(subs_current_tmp) == 0) {
      subs_tmp <- subs_tmp[-1]
      kins_tmp <- kins_tmp[-1]
    }
  } else {
    subs_current_tmp <- subs_current_tmp[-1]
    subs_tmp <- subs_tmp[-1]
    kins_current_tmp <- kins_current_tmp[-1]
    kins_tmp <- kins_tmp[-1]
  }
}