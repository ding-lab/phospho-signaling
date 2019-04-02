# Yige Wu @ March 2019 WashU
# show the number of associated substrate phosphosites per kinase with functional annotation, especially those that haven't been reported before


# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)
source('./cptac2p_analysis/preprocess_files/preprocess_files_shared.R')
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")


# inputs regulatory sites------------------------------------------------------------------
regulatory_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/Mar_04_2019/Regulatory_sites", data.table = F)
regulatory_sites <- regulatory_sites %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1])
regulatory_sites2merge <- regulatory_sites %>%
  select(GENE, SUB_MOD_RSD, ON_FUNCTION, ON_PROCESS, ON_PROT_INTERACT, ON_OTHER_INTERACT) %>%
  unique

# input regression --------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")
regression %>% nrow()
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)
regression %>% nrow()

regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression <- annotate_ks_source(regression = regression)
table(regression$Cancer)


# annotate the regulatory sites -------------------------------------------
regression <- merge(regression, regulatory_sites2merge, by.x = c("SUBSTRATE", "SUB_MOD_RSD"), by.y = c("GENE", "SUB_MOD_RSD"), all.x = T)


# input druggable kinases -------------------------------------------------
drug_genes <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/reference_files/gene_drug_list/Premed_raw_databases/drugBank/drug_list.txt", data.table = F, col.names = "gene")

# get all the terms in ON_PROCESS -----------------------------------------
processes_terms <- regression %>%
  filter(regulated == T) %>%
  select(ON_PROCESS) %>%
  unique() %>%
  unlist()

processes_terms <-  str_split(string = processes_terms, pattern = "; ") %>%
  unlist() %>%
  unique %>%
  sort()

processes_terms_estimate <- regression %>%
  filter(regulated == T) %>%
  select(ON_PROCESS) %>%
  unlist()

processes_terms_estimate <-  str_split(string = processes_terms_estimate, pattern = "; ") %>%
  unlist() %>%
  table()

processes_terms_estimate <- data.frame(ON_PROCESS_term = names(processes_terms_estimate), Freq = as.vector(processes_terms_estimate))
processes_terms_estimate <- processes_terms_estimate %>%
  filter(ON_PROCESS_term != "") %>%
  mutate(ON_PROCESS_term_parent = str_split_fixed(string = ON_PROCESS_term, pattern = ",", n = 2)[,1])
  arrange(Freq)
processes_terms_estimate

# terms2process <- c("endocytosis", "DNA repair", "translation", "cell adhesion", "exocytosis", "autophagy", "cell differentiation", "carcinogenesis", "cell growth", "transcription", "apoptosis", "cell motility")
terms2process <- c("transcription", "apoptosis", "cell growth", "cell differentiation", "cell motility")
# terms2process <- unique(processes_terms_estimate$ON_PROCESS_term_parent)
terms2process
for (ON_PROCESS_term in terms2process) {
  # get table and make figure about ON_PROCESS_tmp -----------------------------------------------
  regression$ON_PROCESS_tmp <- sapply(1:nrow(regression), FUN = function(i, vector_function, ON_PROCESS_term) {
    if (!is.na(vector_function[i]) & grepl(x = vector_function[i], pattern = ON_PROCESS_term)) {
      tmp <- "altered"
      if (!grepl(x = vector_function[i], pattern = paste0(ON_PROCESS_term, ", inhibited")) & grepl(x = vector_function[i], pattern = paste0(ON_PROCESS_term, ", induced"))) {
        tmp <- "induced"
      }
      if (grepl(x = vector_function[i], pattern = paste0(ON_PROCESS_term, ", inhibited")) & !grepl(x = vector_function[i], pattern = paste0(ON_PROCESS_term, ", induced"))) {
        tmp <- "inhibited"
      }
    } else {
      tmp <- NA
    }
    return(tmp)
  }, vector_function = regression$ON_PROCESS, ON_PROCESS_term = ON_PROCESS_term)
  
  kinase_ON_PROCESS_tmp_tab <- regression %>%
    filter(regulated == T) %>%
    filter(SELF == "trans") %>%
    filter(!is.na(ON_PROCESS_tmp)) %>%
    filter(enzyme_type == "kinase") %>%
    mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
    select(GENE, ON_PROCESS_tmp) %>%
    table() %>%
    as.data.frame()
  
  tab2p <- kinase_ON_PROCESS_tmp_tab
  tab2p <- tab2p %>%
    mutate(y = GENE) %>%
    mutate(x = ON_PROCESS_tmp) %>%
    mutate(point_size = Freq) %>%
    mutate(point_fill = ON_PROCESS_tmp) %>%
    filter(Freq > 0) %>%
    filter(x %in% c("induced", "inhibited"))
  tab2p_kinase <- tab2p
  
  ## take out the confusing annotations
  phosphosite2filter <- regression %>%
    filter(regulated == T) %>%
    filter(SELF == "trans") %>%
    filter(!is.na(ON_PROCESS_tmp)) %>%
    filter(GENE %in% tab2p_kinase$GENE) %>%
    filter(ON_PROCESS_tmp %in% c("induced", "inhibited")) %>%
    mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
    select(phosphosite, ON_PROCESS_tmp) %>%
    table()
  if (ncol(phosphosite2filter) == 2) {
    phosphosite2filter <- rownames(phosphosite2filter)[phosphosite2filter[,1] > 0 & phosphosite2filter[,2] > 0]
  } else if (ncol(phosphosite2filter) == 1) {
    next()
    phosphosite2filter <- rownames(phosphosite2filter)[phosphosite2filter[,1] > 0]
  }
  
  kinase_ON_PROCESS_tmp_tab <- regression %>%
    filter(regulated == T) %>%
    filter(SELF == "trans") %>%
    filter(!is.na(ON_PROCESS_tmp)) %>%
    filter(enzyme_type == "kinase") %>%
    mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
    filter(!(phosphosite %in% phosphosite2filter)) %>%
    select(GENE, ON_PROCESS_tmp) %>%
    table() %>%
    as.data.frame()
  
  kinase_order_tab <- regression %>%
    filter(regulated == T) %>%
    filter(SELF == "trans") %>%
    filter(!is.na(ON_PROCESS_tmp)) %>%
    filter(enzyme_type == "kinase") %>%
    filter(ON_PROCESS_tmp %in% c("induced", "inhibited")) %>%
    select(GENE, ON_PROCESS_tmp) %>%
    table() %>%
    as.matrix() 
  if (ncol(kinase_order_tab) == 2) {
    kinase_order_tab <- kinase_order_tab[order(kinase_order_tab[,2], kinase_order_tab[,1], decreasing = T),]
    kinase_order_tab <- rbind(kinase_order_tab[kinase_order_tab[,2] > 0 & kinase_order_tab[,1] == 0,],
                              kinase_order_tab[kinase_order_tab[,2] > 0 & kinase_order_tab[,1] > 0,],
                              kinase_order_tab[kinase_order_tab[,2] == 0 & kinase_order_tab[,1] > 0,])
    kinases_ordered <- rownames(kinase_order_tab)
  } else {
    if (length(kinase_order_tab) <= 1) {
      next()
    }
    kinase_order_tab <- kinase_order_tab[order(kinase_order_tab[,1], decreasing = T),]
    kinases_ordered <- names(kinase_order_tab)
  }
  kinases_ordered <- kinases_ordered[kinases_ordered != ""]
  tab2p$y <- factor(tab2p$y, levels = rownames(kinase_order_tab))
  
  tab4lineplot_y <- tab2p_kinase %>%
    group_by(GENE) %>%
    summarise(Freq_sum = sum(Freq)) %>%
    as.data.frame()
  
  p <- ggplot()
  p <- p + geom_point(data = tab2p, mapping = aes(y = y, x = x, size = point_size, fill = point_fill), shape = 21, color = "white", alpha = 0.5)
  p <- p + geom_text(data = tab2p[tab2p$Freq > 1,], mapping = aes(y = y, x = x, label = Freq))
  p <- p + scale_fill_manual(values = c("altered" = "purple", "induced" = "red", "inhibited" = "blue"))
  p <- p + theme_bw()
  p <- p + theme(axis.title.x = element_blank(), 
                 axis.title.y = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks.x = element_blank(), 
                 axis.ticks.y = element_blank())
  
  fn <- paste0(makeOutDir(resultD = resultD), "num_", ON_PROCESS_term, "_pairs.pdf")
  plot.new()
  pdf(file = fn, width = 2.5, height = 10)
  print(p)
  dev.off()
  
  tab2p <- regression
  tab2p <- tab2p %>%
    filter(regulated == T) %>%
    filter(SELF == "trans") %>%
    filter(!is.na(ON_PROCESS_tmp)) %>%
    filter(GENE %in% tab2p_kinase$GENE) %>%
    filter(ON_PROCESS_tmp %in% c("induced", "inhibited")) %>%
    mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
    filter(!(phosphosite %in% phosphosite2filter)) %>%
    select(GENE, phosphosite, ON_PROCESS_tmp, pair) %>%
    unique
  
  tab2p <- melt(tab2p, id.vars = c("pair", "ON_PROCESS_tmp"))
  tab2p <- tab2p %>%
    # mutate(x = ifelse(grepl(pattern = "_", x = value), "phosphosite", "kinase")) %>%
    arrange(ON_PROCESS_tmp) %>%
    mutate(GENE = str_split_fixed(string = pair, pattern = ":", n = 3)[,1])
  tab2p <- merge(tab2p, tab4lineplot_y, by = c("GENE"), all.x = T)
  
  
  kinase_y <- 0
  phosphosite_y <- 0
  tab2p_new <- NULL
  for (kinase_tmp in kinases_ordered) {
    tab2p_kinase_tmp <- tab2p %>%
      filter(GENE == kinase_tmp) %>%
      filter(variable == "GENE")
    
    tab2p_phosphosite_tmp <- tab2p %>%
      filter(GENE == kinase_tmp) %>%
      filter(variable == "phosphosite")
    
    if (nrow(tab2p_kinase_tmp) == 0 | nrow(tab2p_phosphosite_tmp) == 0) {
      next()
    }
    tab2p_kinase_tmp$y <- kinase_y
    tab2p_new <- rbind(tab2p_new, tab2p_kinase_tmp)
    kinase_y <- kinase_y + 1
    
    y_tmp <- vector(mode = "numeric", length = nrow(tab2p_phosphosite_tmp))
    for (j in 1:nrow(tab2p_phosphosite_tmp)) {
      if (tab2p_phosphosite_tmp$value[j] %in% tab2p_new$value) {
        y_tmp[j] <- tab2p_new$y[tab2p_new$value == tab2p_phosphosite_tmp$value[j]]
      } else {
        y_tmp[j] <- phosphosite_y
        phosphosite_y <- phosphosite_y + 1
      }
    }
    tab2p_phosphosite_tmp$y <- y_tmp
    tab2p_new <- rbind(tab2p_new, tab2p_phosphosite_tmp)
  }
  
  tab2p <- tab2p_new
  tab2p <- tab2p %>%
    mutate(is.direct = (pair %in% regression$pair[regression$is.direct]))
  tab2p$x <- 0
  # tab2p$x[tab2p$variable == "GENE" & tab2p$ON_PROCESS_tmp == "induced"] <- 1
  tab2p$x[tab2p$variable == "GENE" & tab2p$ON_PROCESS_tmp == "induced"] <- -1
  tab2p$x[tab2p$variable == "GENE" & tab2p$ON_PROCESS_tmp == "inhibited"] <- -1
  # tab2p$x[tab2p$variable == "phosphosite" & tab2p$ON_PROCESS_tmp == "induced"] <- 0.2
  tab2p$x[tab2p$variable == "phosphosite" & tab2p$ON_PROCESS_tmp == "induced"] <- -0.2
  tab2p$x[tab2p$variable == "phosphosite" & tab2p$ON_PROCESS_tmp == "inhibited"] <- -0.2

  tmp <- as.vector(tab2p$value)
  tmp[tmp %in% drug_genes$gene] <- paste0("*", tmp[tmp %in% drug_genes$gene])
  tab2p$text <- tmp
  
  p <- ggplot()
  p <- p + geom_label(data = tab2p[tab2p$variable == "GENE" & (tab2p$value %in% tab2p_kinase$GENE[tab2p_kinase$ON_PROCESS_tmp == "induced"])  & !(tab2p$value %in% tab2p_kinase$GENE[tab2p_kinase$ON_PROCESS_tmp == "inhibited"]),], 
                      mapping = aes(x = x, y = y, label = text), nudge_x = -0.2, color = set1[1])
  p <- p + geom_label(data = tab2p[tab2p$variable == "GENE" & !(tab2p$value %in% tab2p_kinase$GENE[tab2p_kinase$ON_PROCESS_tmp == "induced"])  & (tab2p$value %in% tab2p_kinase$GENE[tab2p_kinase$ON_PROCESS_tmp == "inhibited"]),],
                      mapping = aes(x = x, y = y, label = text), nudge_x = -0.2, color = set1[2])
  p <- p + geom_label(data = tab2p[tab2p$variable == "GENE" & (tab2p$value %in% tab2p_kinase$GENE[tab2p_kinase$ON_PROCESS_tmp == "induced"])  & (tab2p$value %in% tab2p_kinase$GENE[tab2p_kinase$ON_PROCESS_tmp == "inhibited"]),],
                      mapping = aes(x = x, y = y, label = text), nudge_x = -0.2, color = "purple")
  p <- p + geom_text(data = tab2p[tab2p$variable == "phosphosite" & tab2p$ON_PROCESS_tmp == "inhibited",], 
                     mapping = aes(x = x, y = y, label = text, color = ON_PROCESS_tmp), nudge_x = 0.2)
  p <- p + geom_text(data = tab2p[tab2p$variable == "phosphosite" & tab2p$ON_PROCESS_tmp == "induced",], 
                     mapping = aes(x = x, y = y, label = text, color = ON_PROCESS_tmp), nudge_x = 0.2)
  p <- p + xlim(c(-1.5, 0.5))
  p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = pair, linetype = is.direct), alpha = 0.6, size = 0.8, color = "black")
  p <- p + scale_color_manual(values = c("FALSE" =  "#FF7F00", "TRUE" = "#6A3D9A", "induced" = set1[1], "inhibited" = set1[2]))
  p <- p + theme_bw() + theme_nogrid()
  p <- p + theme(axis.title.x = element_blank(), 
                 axis.title.y = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks.x = element_blank(), 
                 axis.text.y = element_blank(), 
                 axis.ticks.y = element_blank())
  p
  fn <- paste0(makeOutDir(resultD = resultD), ON_PROCESS_term, "_pairs.pdf")
  plot.new()
  # pdf(file = fn, width = 8, height = 10)
  pdf(file = fn, width = 10, height = 10)
  print(p)
  dev.off()
  
  fn <- paste0(makeOutDir(resultD = resultD), ON_PROCESS_term, "_pairs_long.pdf")
  plot.new()
  # pdf(file = fn, width = 8, height = 10)
  pdf(file = fn, width = 10, height = 15)
  print(p)
  dev.off()
}


regression_dark_kinase <- regression %>%
  filter(regulated == T) %>%
  filter(SELF == "trans") %>%
  filter(GENE %in% kinase_dark_tab$gene[kinase_dark_tab$is_dark_kinase])

regression_dark_kinase %>%
  select(GENE) %>%
  unique()

regression_dark_kinase %>%
  select(GENE) %>%
  unique() %>%
  nrow()

regression %>%
  filter(regulated == T) %>%
  filter(SELF == "trans") %>%
  filter(GENE %in% kinase_dark_tab$gene[kinase_dark_tab$is_drug_kinase]) %>%
  select(GENE) %>%
  unique() %>%
  nrow()
