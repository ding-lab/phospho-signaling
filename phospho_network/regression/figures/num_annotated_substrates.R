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

# regression <- change_regression_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression <- annotate_ks_source(regression = regression)
table(regression$Cancer)


# annotate the regulatory sites -------------------------------------------
regression <- merge(regression, regulatory_sites2merge, by.x = c("SUBSTRATE", "SUB_MOD_RSD"), by.y = c("GENE", "SUB_MOD_RSD"), all.x = T)


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
  arrange(Freq)
processes_terms_estimate

# get table and make figure about apoptosis -----------------------------------------------
regression$apoptosis <- sapply(1:nrow(regression), FUN = function(i, vector_function) {
  if (!is.na(vector_function[i]) & grepl(x = vector_function[i], pattern = "apoptosis")) {
    tmp <- "altered"
    if (!grepl(x = vector_function[i], pattern = "apoptosis, inhibited") & grepl(x = vector_function[i], pattern = "apoptosis, induced")) {
      tmp <- "induced"
    }
    if (grepl(x = vector_function[i], pattern = "apoptosis, inhibited") & !grepl(x = vector_function[i], pattern = "apoptosis, induced")) {
      tmp <- "inhibited"
    }
  } else {
    tmp <- NA
  }
  return(tmp)
}, vector_function = regression$ON_PROCESS)
kinase_apoptosis_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(apoptosis)) %>%
  filter(enzyme_type == "kinase") %>%
  select(GENE, apoptosis) %>%
  table() %>%
  as.data.frame()

kinase_order_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(apoptosis)) %>%
  filter(enzyme_type == "kinase") %>%
  filter(apoptosis %in% c("induced", "inhibited")) %>%
  select(GENE, apoptosis) %>%
  table() %>%
  as.matrix() 
kinase_order_tab <- kinase_order_tab[order(kinase_order_tab[,2], kinase_order_tab[,1], decreasing = T),]

tab2p <- kinase_apoptosis_tab
tab2p <- tab2p %>%
  mutate(y = GENE) %>%
  mutate(x = apoptosis) %>%
  mutate(point_size = Freq) %>%
  mutate(point_fill = apoptosis) %>%
  filter(Freq > 0) %>%
  filter(x %in% c("induced", "inhibited"))
tab2p$y <- factor(tab2p$y, levels = rownames(kinase_order_tab))
tab2p_kinase <- tab2p

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

fn <- paste0(makeOutDir(resultD = resultD), "num_apoptosis_pairs.pdf")
pdf(file = fn, width = 2.5, height = 10)
p
dev.off()

tab2p <- regression
tab2p <- tab2p %>%
  filter(regulated == T) %>%
  filter(!is.na(apoptosis)) %>%
  filter(GENE %in% tab2p_kinase$GENE) %>%
  filter(apoptosis %in% c("induced", "inhibited")) %>%
  mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
  select(GENE, phosphosite, apoptosis, pair) %>%
  unique

tab2p <- melt(tab2p, id.vars = c("pair", "apoptosis"))
tab2p <- tab2p %>%
  mutate(x = ifelse(grepl(pattern = "_", x = value), "phosphosite", "kinase")) %>%
  arrange(apoptosis) %>%
  mutate(GENE = str_split_fixed(string = pair, pattern = ":", n = 3)[,1])
tab2p <- merge(tab2p, tab4lineplot_y, by = c("GENE"), all.x = T)


kinase_y <- 0
phosphosite_y <- 0
tab2p_new <- NULL
for (kinase_tmp in rownames(kinase_order_tab)) {
  tab2p_kinase_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "GENE")
  
  tab2p_phosphosite_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "phosphosite")
  
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

p <- ggplot()
p <- p + geom_text(data = tab2p[tab2p$x == "kinase",], mapping = aes(x = x, y = y, label = value), nudge_x = -0.2)
p <- p + geom_text(data = tab2p[tab2p$x == "phosphosite",], mapping = aes(x = x, y = y, label = value, color = apoptosis), 
                   nudge_x = 0.3)
p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = pair, color = apoptosis, linetype = is.direct), alpha = 0.5)
p <- p + scale_color_manual(values = c("FALSE" =  "orange", "TRUE" = "grey50", "induced" = set1[1], "inhibited" = set1[2]))
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank(), 
               axis.title.y = element_blank(), 
               axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(), 
               axis.text.y = element_blank(), 
               axis.ticks.y = element_blank())
fn <- paste0(makeOutDir(resultD = resultD), "apoptosis_pairs.pdf")
pdf(file = fn, width = 6, height = 10)
p
dev.off()

# get table and make figure about transcription -----------------------------------------------
regression$transcription <- sapply(1:nrow(regression), FUN = function(i, vector_function) {
  if (!is.na(vector_function[i]) & grepl(x = vector_function[i], pattern = "transcription")) {
    tmp <- "altered"
    if (!grepl(x = vector_function[i], pattern = "transcription, inhibited") & grepl(x = vector_function[i], pattern = "transcription, induced")) {
      tmp <- "induced"
    }
    if (grepl(x = vector_function[i], pattern = "transcription, inhibited") & !grepl(x = vector_function[i], pattern = "transcription, induced")) {
      tmp <- "inhibited"
    }
  } else {
    tmp <- NA
  }
  return(tmp)
}, vector_function = regression$ON_PROCESS)
kinase_transcription_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(transcription)) %>%
  filter(enzyme_type == "kinase") %>%
  select(GENE, transcription) %>%
  table() %>%
  as.data.frame()
kinase_transcription_tab


kinase_order_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(transcription)) %>%
  filter(enzyme_type == "kinase") %>%
  filter(transcription %in% c("induced", "inhibited")) %>%
  select(GENE, transcription) %>%
  table() %>%
  as.matrix() 
kinase_order_tab <- kinase_order_tab[order(kinase_order_tab[,2], kinase_order_tab[,1], decreasing = T),]

tab2p <- kinase_transcription_tab
tab2p <- tab2p %>%
  mutate(y = GENE) %>%
  mutate(x = transcription) %>%
  mutate(point_size = Freq) %>%
  mutate(point_fill = transcription) %>%
  filter(Freq > 0) %>%
  filter(x %in% c("induced", "inhibited"))
tab2p$y <- factor(tab2p$y, levels = rownames(kinase_order_tab))
tab2p_kinase <- tab2p

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

fn <- paste0(makeOutDir(resultD = resultD), "num_transcription_pairs.pdf")
pdf(file = fn, width = 2.5, height = 10)
p
dev.off()

tab2p <- regression
tab2p <- tab2p %>%
  filter(regulated == T) %>%
  filter(!is.na(transcription)) %>%
  filter(GENE %in% tab2p_kinase$GENE) %>%
  filter(transcription %in% c("induced", "inhibited")) %>%
  mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
  select(GENE, phosphosite, transcription, pair) %>%
  unique

tab2p <- melt(tab2p, id.vars = c("pair", "transcription"))
tab2p <- tab2p %>%
  mutate(x = ifelse(grepl(pattern = "_", x = value), "phosphosite", "kinase")) %>%
  arrange(transcription) %>%
  mutate(GENE = str_split_fixed(string = pair, pattern = ":", n = 3)[,1])
tab2p <- merge(tab2p, tab4lineplot_y, by = c("GENE"), all.x = T)


kinase_y <- 0
phosphosite_y <- 0
tab2p_new <- NULL
for (kinase_tmp in rownames(kinase_order_tab)) {
  tab2p_kinase_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "GENE")
  
  tab2p_phosphosite_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "phosphosite")
  
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

p <- ggplot()
p <- p + geom_text(data = tab2p[tab2p$x == "kinase",], mapping = aes(x = x, y = y, label = value), nudge_x = -0.2)
p <- p + geom_text(data = tab2p[tab2p$x == "phosphosite",], mapping = aes(x = x, y = y, label = value, color = transcription), 
                   nudge_x = 0.3)
p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = pair, color = transcription, linetype = is.direct), alpha = 0.5)
p <- p + scale_color_manual(values = c("FALSE" =  "orange", "TRUE" = "grey50", "induced" = set1[1], "inhibited" = set1[2]))
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank(), 
               axis.title.y = element_blank(), 
               axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(), 
               axis.text.y = element_blank(), 
               axis.ticks.y = element_blank())
fn <- paste0(makeOutDir(resultD = resultD), "transcription_pairs.pdf")
pdf(file = fn, width = 6, height = 10)
p
dev.off()




# get table and make figure about cell motility -----------------------------------------------
regression$cell_motility <- sapply(1:nrow(regression), FUN = function(i, vector_function) {
  if (!is.na(vector_function[i]) & grepl(x = vector_function[i], pattern = "cell motility")) {
    tmp <- "altered"
    if (!grepl(x = vector_function[i], pattern = "cell motility, inhibited") & grepl(x = vector_function[i], pattern = "cell motility, induced")) {
      tmp <- "induced"
    }
    if (grepl(x = vector_function[i], pattern = "cell motility, inhibited") & !grepl(x = vector_function[i], pattern = "cell motility, induced")) {
      tmp <- "inhibited"
    }
  } else {
    tmp <- NA
  }
  return(tmp)
}, vector_function = regression$ON_PROCESS)
kinase_cell_motility_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(cell_motility)) %>%
  filter(enzyme_type == "kinase") %>%
  select(GENE, cell_motility) %>%
  table() %>%
  as.data.frame()
kinase_cell_motility_tab

kinase_order_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(cell_motility)) %>%
  filter(enzyme_type == "kinase") %>%
  filter(cell_motility %in% c("induced", "inhibited")) %>%
  select(GENE, cell_motility) %>%
  table() %>%
  as.matrix()
kinase_order_tab <- kinase_order_tab[order(kinase_order_tab[,2], kinase_order_tab[,1], decreasing = T),]

tab2p <- kinase_cell_motility_tab
tab2p <- tab2p %>%
  mutate(y = GENE) %>%
  mutate(x = cell_motility) %>%
  mutate(point_size = Freq) %>%
  mutate(point_fill = cell_motility) %>%
  filter(Freq > 0) %>%
  filter(x %in% c("induced", "inhibited"))
tab2p$y <- factor(tab2p$y, levels = rownames(kinase_order_tab))
tab2p_kinase <- tab2p

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

fn <- paste0(makeOutDir(resultD = resultD), "num_cell_motility_pairs.pdf")
pdf(file = fn, width = 2.5, height = 10)
p
dev.off()

tab2p <- regression
tab2p <- tab2p %>%
  filter(regulated == T) %>%
  filter(!is.na(cell_motility)) %>%
  filter(GENE %in% tab2p_kinase$GENE) %>%
  filter(cell_motility %in% c("induced", "inhibited")) %>%
  mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
  select(GENE, phosphosite, cell_motility, pair) %>%
  unique

tab2p <- melt(tab2p, id.vars = c("pair", "cell_motility"))
tab2p <- tab2p %>%
  mutate(x = ifelse(grepl(pattern = "_", x = value), "phosphosite", "kinase")) %>%
  arrange(cell_motility) %>%
  mutate(GENE = str_split_fixed(string = pair, pattern = ":", n = 3)[,1])
tab2p <- merge(tab2p, tab4lineplot_y, by = c("GENE"), all.x = T)


kinase_y <- 0
phosphosite_y <- 0
tab2p_new <- NULL
for (kinase_tmp in rownames(kinase_order_tab)) {
  tab2p_kinase_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "GENE")
  
  tab2p_phosphosite_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "phosphosite")
  
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

p <- ggplot()
p <- p + geom_text(data = tab2p[tab2p$x == "kinase",], mapping = aes(x = x, y = y, label = value), nudge_x = -0.2)
p <- p + geom_text(data = tab2p[tab2p$x == "phosphosite",], mapping = aes(x = x, y = y, label = value, color = cell_motility), 
                   nudge_x = 0.3)
p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = pair, color = cell_motility, linetype = is.direct), alpha = 0.5)
p <- p + scale_color_manual(values = c("FALSE" =  "orange", "TRUE" = "grey50", "induced" = set1[1], "inhibited" = set1[2]))
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank(), 
               axis.title.y = element_blank(), 
               axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(), 
               axis.text.y = element_blank(), 
               axis.ticks.y = element_blank())
fn <- paste0(makeOutDir(resultD = resultD), "cell_motility_pairs.pdf")
pdf(file = fn, width = 6, height = 10)
p
dev.off()















# get table and make figure about cell growth -----------------------------------------------
regression$cell_growth <- sapply(1:nrow(regression), FUN = function(i, vector_function) {
  if (!is.na(vector_function[i]) & grepl(x = vector_function[i], pattern = "cell growth")) {
    tmp <- "altered"
    if (!grepl(x = vector_function[i], pattern = "cell growth, inhibited") & grepl(x = vector_function[i], pattern = "cell growth, induced")) {
      tmp <- "induced"
    }
    if (grepl(x = vector_function[i], pattern = "cell growth, inhibited") & !grepl(x = vector_function[i], pattern = "cell growth, induced")) {
      tmp <- "inhibited"
    }
  } else {
    tmp <- NA
  }
  return(tmp)
}, vector_function = regression$ON_PROCESS)


kinase_cell_growth_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(cell_growth)) %>%
  filter(enzyme_type == "kinase") %>%
  select(GENE, cell_growth) %>%
  table() %>%
  as.data.frame()
kinase_cell_growth_tab

kinase_order_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(cell_growth)) %>%
  filter(enzyme_type == "kinase") %>%
  filter(cell_growth %in% c("induced", "inhibited")) %>%
  select(GENE, cell_growth) %>%
  table() %>%
  as.matrix()
kinase_order_tab <- kinase_order_tab[order(kinase_order_tab[,2], kinase_order_tab[,1], decreasing = T),]

tab2p <- kinase_cell_growth_tab
tab2p <- tab2p %>%
  mutate(y = GENE) %>%
  mutate(x = cell_growth) %>%
  mutate(point_size = Freq) %>%
  mutate(point_fill = cell_growth) %>%
  filter(Freq > 0) %>%
  filter(x %in% c("induced", "inhibited"))
tab2p$y <- factor(tab2p$y, levels = rownames(kinase_order_tab))
tab2p_kinase <- tab2p

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

fn <- paste0(makeOutDir(resultD = resultD), "num_cell_growth_pairs.pdf")
pdf(file = fn, width = 2.5, height = 10)
p
dev.off()

tab2p <- regression
tab2p <- tab2p %>%
  filter(regulated == T) %>%
  filter(!is.na(cell_growth)) %>%
  filter(GENE %in% tab2p_kinase$GENE) %>%
  filter(cell_growth %in% c("induced", "inhibited")) %>%
  mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
  select(GENE, phosphosite, cell_growth, pair) %>%
  unique

tab2p <- melt(tab2p, id.vars = c("pair", "cell_growth"))
tab2p <- tab2p %>%
  mutate(x = ifelse(grepl(pattern = "_", x = value), "phosphosite", "kinase")) %>%
  arrange(cell_growth) %>%
  mutate(GENE = str_split_fixed(string = pair, pattern = ":", n = 3)[,1])
tab2p <- merge(tab2p, tab4lineplot_y, by = c("GENE"), all.x = T)


kinase_y <- 0
phosphosite_y <- 0
tab2p_new <- NULL
for (kinase_tmp in rownames(kinase_order_tab)) {
  tab2p_kinase_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "GENE")
  
  tab2p_phosphosite_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "phosphosite")
  
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

p <- ggplot()
p <- p + geom_text(data = tab2p[tab2p$x == "kinase",], mapping = aes(x = x, y = y, label = value), nudge_x = -0.2)
p <- p + geom_text(data = tab2p[tab2p$x == "phosphosite",], mapping = aes(x = x, y = y, label = value, color = cell_growth), 
                   nudge_x = 0.3)
p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = pair, color = cell_growth), alpha = 0.5)
p <- p + scale_color_manual(values = c("FALSE" =  "orange", "TRUE" = "grey50", "induced" = set1[1], "inhibited" = set1[2]))
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank(), 
               axis.title.y = element_blank(), 
               axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(), 
               axis.text.y = element_blank(), 
               axis.ticks.y = element_blank())
fn <- paste0(makeOutDir(resultD = resultD), "cell_growth_pairs.pdf")
plot.new()
pdf(file = fn, width = 6, height = 10)
p
dev.off()





# get table and make figure about carcinogenesis -----------------------------------------------
regression$carcinogenesis <- sapply(1:nrow(regression), FUN = function(i, vector_function) {
  if (!is.na(vector_function[i]) & grepl(x = vector_function[i], pattern = "carcinogenesis")) {
    tmp <- "altered"
    if (!grepl(x = vector_function[i], pattern = "carcinogenesis, inhibited") & grepl(x = vector_function[i], pattern = "carcinogenesis, induced")) {
      tmp <- "induced"
    }
    if (grepl(x = vector_function[i], pattern = "carcinogenesis, inhibited") & !grepl(x = vector_function[i], pattern = "carcinogenesis, induced")) {
      tmp <- "inhibited"
    }
  } else {
    tmp <- NA
  }
  return(tmp)
}, vector_function = regression$ON_PROCESS)
kinase_carcinogenesis_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(carcinogenesis)) %>%
  filter(enzyme_type == "kinase") %>%
  select(GENE, carcinogenesis) %>%
  table() %>%
  as.data.frame()
kinase_carcinogenesis_tab

kinase_order_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(carcinogenesis)) %>%
  filter(enzyme_type == "kinase") %>%
  filter(carcinogenesis %in% c("induced", "inhibited")) %>%
  select(GENE, carcinogenesis) %>%
  table() %>%
  as.matrix()
kinase_order_tab <- kinase_order_tab[order(kinase_order_tab[,2], kinase_order_tab[,1], decreasing = T),]

tab2p <- kinase_carcinogenesis_tab
tab2p <- tab2p %>%
  mutate(y = GENE) %>%
  mutate(x = carcinogenesis) %>%
  mutate(point_size = Freq) %>%
  mutate(point_fill = carcinogenesis) %>%
  filter(Freq > 0) %>%
  filter(x %in% c("induced", "inhibited"))
tab2p$y <- factor(tab2p$y, levels = rownames(kinase_order_tab))
tab2p_kinase <- tab2p

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

fn <- paste0(makeOutDir(resultD = resultD), "num_carcinogenesis_pairs.pdf")
pdf(file = fn, width = 2.5, height = 10)
p
dev.off()

tab2p <- regression
tab2p <- tab2p %>%
  filter(regulated == T) %>%
  filter(!is.na(carcinogenesis)) %>%
  filter(GENE %in% tab2p_kinase$GENE) %>%
  filter(carcinogenesis %in% c("induced", "inhibited")) %>%
  mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
  select(GENE, phosphosite, carcinogenesis, pair) %>%
  unique

tab2p <- melt(tab2p, id.vars = c("pair", "carcinogenesis"))
tab2p <- tab2p %>%
  mutate(x = ifelse(grepl(pattern = "_", x = value), "phosphosite", "kinase")) %>%
  arrange(carcinogenesis) %>%
  mutate(GENE = str_split_fixed(string = pair, pattern = ":", n = 3)[,1])
tab2p <- merge(tab2p, tab4lineplot_y, by = c("GENE"), all.x = T)


kinase_y <- 0
phosphosite_y <- 0
tab2p_new <- NULL
for (kinase_tmp in rownames(kinase_order_tab)) {
  tab2p_kinase_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "GENE")
  
  tab2p_phosphosite_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "phosphosite")
  
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

p <- ggplot()
p <- p + geom_text(data = tab2p[tab2p$x == "kinase",], mapping = aes(x = x, y = y, label = value), nudge_x = -0.2)
p <- p + geom_text(data = tab2p[tab2p$x == "phosphosite",], mapping = aes(x = x, y = y, label = value, color = carcinogenesis), 
                   nudge_x = 0.3)
p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = pair, color = carcinogenesis), alpha = 0.5)
p <- p + scale_color_manual(values = c("FALSE" =  "orange", "TRUE" = "grey50", "induced" = set1[1], "inhibited" = set1[2]))
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank(), 
               axis.title.y = element_blank(), 
               axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(), 
               axis.text.y = element_blank(), 
               axis.ticks.y = element_blank())
fn <- paste0(makeOutDir(resultD = resultD), "carcinogenesis_pairs.pdf")
plot.new()
pdf(file = fn, width = 6, height = 10)
p
dev.off()







# get table and make figure about cell differentiation -----------------------------------------------
regression$cell_differentiation <- sapply(1:nrow(regression), FUN = function(i, vector_function) {
  if (!is.na(vector_function[i]) & grepl(x = vector_function[i], pattern = "cell differentiation")) {
    tmp <- "altered"
    if (!grepl(x = vector_function[i], pattern = "cell differentiation, inhibited") & grepl(x = vector_function[i], pattern = "cell differentiation, induced")) {
      tmp <- "induced"
    }
    if (grepl(x = vector_function[i], pattern = "cell differentiation, inhibited") & !grepl(x = vector_function[i], pattern = "cell differentiation, induced")) {
      tmp <- "inhibited"
    }
  } else {
    tmp <- NA
  }
  return(tmp)
}, vector_function = regression$ON_PROCESS)
kinase_cell_differentiation_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(cell_differentiation)) %>%
  filter(enzyme_type == "kinase") %>%
  select(GENE, cell_differentiation) %>%
  table() %>%
  as.data.frame()
kinase_cell_differentiation_tab

kinase_order_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(cell_differentiation)) %>%
  filter(enzyme_type == "kinase") %>%
  filter(cell_differentiation %in% c("induced", "inhibited")) %>%
  select(GENE, cell_differentiation) %>%
  table() %>%
  as.matrix()
kinase_order_tab <- kinase_order_tab[order(kinase_order_tab[,2], kinase_order_tab[,1], decreasing = T),]

tab2p <- kinase_cell_differentiation_tab
tab2p <- tab2p %>%
  mutate(y = GENE) %>%
  mutate(x = cell_differentiation) %>%
  mutate(point_size = Freq) %>%
  mutate(point_fill = cell_differentiation) %>%
  filter(Freq > 0) %>%
  filter(x %in% c("induced", "inhibited"))
tab2p$y <- factor(tab2p$y, levels = rownames(kinase_order_tab))
tab2p_kinase <- tab2p

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

fn <- paste0(makeOutDir(resultD = resultD), "num_cell_differentiation_pairs.pdf")
pdf(file = fn, width = 2.5, height = 10)
p
dev.off()

tab2p <- regression
tab2p <- tab2p %>%
  filter(regulated == T) %>%
  filter(!is.na(cell_differentiation)) %>%
  filter(GENE %in% tab2p_kinase$GENE) %>%
  filter(cell_differentiation %in% c("induced", "inhibited")) %>%
  mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
  select(GENE, phosphosite, cell_differentiation, pair) %>%
  unique

tab2p <- melt(tab2p, id.vars = c("pair", "cell_differentiation"))
tab2p <- tab2p %>%
  mutate(x = ifelse(grepl(pattern = "_", x = value), "phosphosite", "kinase")) %>%
  arrange(cell_differentiation) %>%
  mutate(GENE = str_split_fixed(string = pair, pattern = ":", n = 3)[,1])
tab2p <- merge(tab2p, tab4lineplot_y, by = c("GENE"), all.x = T)


kinase_y <- 0
phosphosite_y <- 0
tab2p_new <- NULL
for (kinase_tmp in rownames(kinase_order_tab)) {
  tab2p_kinase_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "GENE")
  
  tab2p_phosphosite_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "phosphosite")
  
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

p <- ggplot()
p <- p + geom_text(data = tab2p[tab2p$x == "kinase",], mapping = aes(x = x, y = y, label = value), nudge_x = -0.2)
p <- p + geom_text(data = tab2p[tab2p$x == "phosphosite",], mapping = aes(x = x, y = y, label = value, color = cell_differentiation), 
                   nudge_x = 0.3)
p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = pair, color = cell_differentiation), alpha = 0.5)
p <- p + scale_color_manual(values = c("FALSE" =  "orange", "TRUE" = "grey50", "induced" = set1[1], "inhibited" = set1[2]))
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank(), 
               axis.title.y = element_blank(), 
               axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(), 
               axis.text.y = element_blank(), 
               axis.ticks.y = element_blank())
fn <- paste0(makeOutDir(resultD = resultD), "cell_differentiation_pairs.pdf")
plot.new()
pdf(file = fn, width = 6, height = 10)
p
dev.off()







# get table and make figure about cell adhesion -----------------------------------------------
regression$cell_adhesion <- sapply(1:nrow(regression), FUN = function(i, vector_function) {
  if (!is.na(vector_function[i]) & grepl(x = vector_function[i], pattern = "cell adhesion")) {
    tmp <- "altered"
    if (!grepl(x = vector_function[i], pattern = "cell adhesion, inhibited") & grepl(x = vector_function[i], pattern = "cell adhesion, induced")) {
      tmp <- "induced"
    }
    if (grepl(x = vector_function[i], pattern = "cell adhesion, inhibited") & !grepl(x = vector_function[i], pattern = "cell adhesion, induced")) {
      tmp <- "inhibited"
    }
  } else {
    tmp <- NA
  }
  return(tmp)
}, vector_function = regression$ON_PROCESS)
kinase_cell_adhesion_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(cell_adhesion)) %>%
  filter(enzyme_type == "kinase") %>%
  select(GENE, cell_adhesion) %>%
  table() %>%
  as.data.frame()
kinase_cell_adhesion_tab

kinase_order_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(cell_adhesion)) %>%
  filter(enzyme_type == "kinase") %>%
  filter(cell_adhesion %in% c("induced", "inhibited")) %>%
  select(GENE, cell_adhesion) %>%
  table() %>%
  as.matrix()
kinase_order_tab <- kinase_order_tab[order(kinase_order_tab[,2], kinase_order_tab[,1], decreasing = T),]

tab2p <- kinase_cell_adhesion_tab
tab2p <- tab2p %>%
  mutate(y = GENE) %>%
  mutate(x = cell_adhesion) %>%
  mutate(point_size = Freq) %>%
  mutate(point_fill = cell_adhesion) %>%
  filter(Freq > 0) %>%
  filter(x %in% c("induced", "inhibited"))
tab2p$y <- factor(tab2p$y, levels = rownames(kinase_order_tab))
tab2p_kinase <- tab2p

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

fn <- paste0(makeOutDir(resultD = resultD), "num_cell_adhesion_pairs.pdf")
pdf(file = fn, width = 2.5, height = 10)
p
dev.off()

tab2p <- regression
tab2p <- tab2p %>%
  filter(regulated == T) %>%
  filter(!is.na(cell_adhesion)) %>%
  filter(GENE %in% tab2p_kinase$GENE) %>%
  filter(cell_adhesion %in% c("induced", "inhibited")) %>%
  mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
  select(GENE, phosphosite, cell_adhesion, pair) %>%
  unique

tab2p <- melt(tab2p, id.vars = c("pair", "cell_adhesion"))
tab2p <- tab2p %>%
  mutate(x = ifelse(grepl(pattern = "_", x = value), "phosphosite", "kinase")) %>%
  arrange(cell_adhesion) %>%
  mutate(GENE = str_split_fixed(string = pair, pattern = ":", n = 3)[,1])
tab2p <- merge(tab2p, tab4lineplot_y, by = c("GENE"), all.x = T)


kinase_y <- 0
phosphosite_y <- 0
tab2p_new <- NULL
for (kinase_tmp in rownames(kinase_order_tab)) {
  tab2p_kinase_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "GENE")
  
  tab2p_phosphosite_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "phosphosite")
  
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

p <- ggplot()
p <- p + geom_text(data = tab2p[tab2p$x == "kinase",], mapping = aes(x = x, y = y, label = value), nudge_x = -0.2)
p <- p + geom_text(data = tab2p[tab2p$x == "phosphosite",], mapping = aes(x = x, y = y, label = value, color = cell_adhesion), 
                   nudge_x = 0.3)
p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = pair, color = cell_adhesion), alpha = 0.5)
p <- p + scale_color_manual(values = c("FALSE" =  "orange", "TRUE" = "grey50", "induced" = set1[1], "inhibited" = set1[2]))
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank(), 
               axis.title.y = element_blank(), 
               axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(), 
               axis.text.y = element_blank(), 
               axis.ticks.y = element_blank())
fn <- paste0(makeOutDir(resultD = resultD), "cell_adhesion_pairs.pdf")
plot.new()
pdf(file = fn, width = 6, height = 10)
p
dev.off()









# get table and make figure about autophagy -----------------------------------------------
regression$autophagy <- sapply(1:nrow(regression), FUN = function(i, vector_function) {
  if (!is.na(vector_function[i]) & grepl(x = vector_function[i], pattern = "autophagy")) {
    tmp <- "altered"
    if (!grepl(x = vector_function[i], pattern = "autophagy, inhibited") & grepl(x = vector_function[i], pattern = "autophagy, induced")) {
      tmp <- "induced"
    }
    if (grepl(x = vector_function[i], pattern = "autophagy, inhibited") & !grepl(x = vector_function[i], pattern = "autophagy, induced")) {
      tmp <- "inhibited"
    }
  } else {
    tmp <- NA
  }
  return(tmp)
}, vector_function = regression$ON_PROCESS)
kinase_autophagy_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(autophagy)) %>%
  filter(enzyme_type == "kinase") %>%
  select(GENE, autophagy) %>%
  table() %>%
  as.data.frame()
kinase_autophagy_tab

kinase_order_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(autophagy)) %>%
  filter(enzyme_type == "kinase") %>%
  filter(autophagy %in% c("induced", "inhibited")) %>%
  select(GENE, autophagy) %>%
  table() %>%
  as.matrix()
kinase_order_tab <- kinase_order_tab[order(kinase_order_tab[,2], kinase_order_tab[,1], decreasing = T),]

tab2p <- kinase_autophagy_tab
tab2p <- tab2p %>%
  mutate(y = GENE) %>%
  mutate(x = autophagy) %>%
  mutate(point_size = Freq) %>%
  mutate(point_fill = autophagy) %>%
  filter(Freq > 0) %>%
  filter(x %in% c("induced", "inhibited"))
tab2p$y <- factor(tab2p$y, levels = rownames(kinase_order_tab))
tab2p_kinase <- tab2p

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

fn <- paste0(makeOutDir(resultD = resultD), "num_autophagy_pairs.pdf")
pdf(file = fn, width = 2.5, height = 10)
p
dev.off()

tab2p <- regression
tab2p <- tab2p %>%
  filter(regulated == T) %>%
  filter(!is.na(autophagy)) %>%
  filter(GENE %in% tab2p_kinase$GENE) %>%
  filter(autophagy %in% c("induced", "inhibited")) %>%
  mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
  select(GENE, phosphosite, autophagy, pair) %>%
  unique

tab2p <- melt(tab2p, id.vars = c("pair", "autophagy"))
tab2p <- tab2p %>%
  mutate(x = ifelse(grepl(pattern = "_", x = value), "phosphosite", "kinase")) %>%
  arrange(autophagy) %>%
  mutate(GENE = str_split_fixed(string = pair, pattern = ":", n = 3)[,1])
tab2p <- merge(tab2p, tab4lineplot_y, by = c("GENE"), all.x = T)


kinase_y <- 0
phosphosite_y <- 0
tab2p_new <- NULL
for (kinase_tmp in rownames(kinase_order_tab)) {
  tab2p_kinase_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "GENE")
  
  tab2p_phosphosite_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "phosphosite")
  
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

p <- ggplot()
p <- p + geom_text(data = tab2p[tab2p$x == "kinase",], mapping = aes(x = x, y = y, label = value), nudge_x = -0.2)
p <- p + geom_text(data = tab2p[tab2p$x == "phosphosite",], mapping = aes(x = x, y = y, label = value, color = autophagy), 
                   nudge_x = 0.3)
p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = pair, color = autophagy), alpha = 0.5)
p <- p + scale_color_manual(values = c("FALSE" =  "orange", "TRUE" = "grey50", "induced" = set1[1], "inhibited" = set1[2]))
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank(), 
               axis.title.y = element_blank(), 
               axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(), 
               axis.text.y = element_blank(), 
               axis.ticks.y = element_blank())
fn <- paste0(makeOutDir(resultD = resultD), "autophagy_pairs.pdf")
plot.new()
pdf(file = fn, width = 6, height = 10)
p
dev.off()









# get table and make figure about exocytosis -----------------------------------------------
regression$exocytosis <- sapply(1:nrow(regression), FUN = function(i, vector_function) {
  if (!is.na(vector_function[i]) & grepl(x = vector_function[i], pattern = "exocytosis")) {
    tmp <- "altered"
    if (!grepl(x = vector_function[i], pattern = "exocytosis, inhibited") & grepl(x = vector_function[i], pattern = "exocytosis, induced")) {
      tmp <- "induced"
    }
    if (grepl(x = vector_function[i], pattern = "exocytosis, inhibited") & !grepl(x = vector_function[i], pattern = "exocytosis, induced")) {
      tmp <- "inhibited"
    }
  } else {
    tmp <- NA
  }
  return(tmp)
}, vector_function = regression$ON_PROCESS)
kinase_exocytosis_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(exocytosis)) %>%
  filter(enzyme_type == "kinase") %>%
  select(GENE, exocytosis) %>%
  table() %>%
  as.data.frame()
kinase_exocytosis_tab

kinase_order_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(exocytosis)) %>%
  filter(enzyme_type == "kinase") %>%
  filter(exocytosis %in% c("induced", "inhibited")) %>%
  select(GENE, exocytosis) %>%
  table() %>%
  as.matrix()
kinase_order_tab <- kinase_order_tab[order(kinase_order_tab[,2], kinase_order_tab[,1], decreasing = T),]

tab2p <- kinase_exocytosis_tab
tab2p <- tab2p %>%
  mutate(y = GENE) %>%
  mutate(x = exocytosis) %>%
  mutate(point_size = Freq) %>%
  mutate(point_fill = exocytosis) %>%
  filter(Freq > 0) %>%
  filter(x %in% c("induced", "inhibited"))
tab2p$y <- factor(tab2p$y, levels = rownames(kinase_order_tab))
tab2p_kinase <- tab2p

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

fn <- paste0(makeOutDir(resultD = resultD), "num_exocytosis_pairs.pdf")
pdf(file = fn, width = 2.5, height = 10)
p
dev.off()

tab2p <- regression
tab2p <- tab2p %>%
  filter(regulated == T) %>%
  filter(!is.na(exocytosis)) %>%
  filter(GENE %in% tab2p_kinase$GENE) %>%
  filter(exocytosis %in% c("induced", "inhibited")) %>%
  mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
  select(GENE, phosphosite, exocytosis, pair) %>%
  unique

tab2p <- melt(tab2p, id.vars = c("pair", "exocytosis"))
tab2p <- tab2p %>%
  mutate(x = ifelse(grepl(pattern = "_", x = value), "phosphosite", "kinase")) %>%
  arrange(exocytosis) %>%
  mutate(GENE = str_split_fixed(string = pair, pattern = ":", n = 3)[,1])
tab2p <- merge(tab2p, tab4lineplot_y, by = c("GENE"), all.x = T)


kinase_y <- 0
phosphosite_y <- 0
tab2p_new <- NULL
for (kinase_tmp in rownames(kinase_order_tab)) {
  tab2p_kinase_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "GENE")
  
  tab2p_phosphosite_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "phosphosite")
  
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

p <- ggplot()
p <- p + geom_text(data = tab2p[tab2p$x == "kinase",], mapping = aes(x = x, y = y, label = value), nudge_x = -0.2)
p <- p + geom_text(data = tab2p[tab2p$x == "phosphosite",], mapping = aes(x = x, y = y, label = value, color = exocytosis), 
                   nudge_x = 0.3)
p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = pair, color = exocytosis), alpha = 0.5)
p <- p + scale_color_manual(values = c("FALSE" =  "orange", "TRUE" = "grey50", "induced" = set1[1], "inhibited" = set1[2]))
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank(), 
               axis.title.y = element_blank(), 
               axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(), 
               axis.text.y = element_blank(), 
               axis.ticks.y = element_blank())
fn <- paste0(makeOutDir(resultD = resultD), "exocytosis_pairs.pdf")
plot.new()
pdf(file = fn, width = 6, height = 10)
p
dev.off()











# get table and make figure about translation -----------------------------------------------
regression$translation <- sapply(1:nrow(regression), FUN = function(i, vector_function) {
  if (!is.na(vector_function[i]) & grepl(x = vector_function[i], pattern = "translation")) {
    tmp <- "altered"
    if (!grepl(x = vector_function[i], pattern = "translation, inhibited") & grepl(x = vector_function[i], pattern = "translation, induced")) {
      tmp <- "induced"
    }
    if (grepl(x = vector_function[i], pattern = "translation, inhibited") & !grepl(x = vector_function[i], pattern = "translation, induced")) {
      tmp <- "inhibited"
    }
  } else {
    tmp <- NA
  }
  return(tmp)
}, vector_function = regression$ON_PROCESS)
kinase_translation_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(translation)) %>%
  filter(enzyme_type == "kinase") %>%
  select(GENE, translation) %>%
  table() %>%
  as.data.frame()
kinase_translation_tab

kinase_order_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(translation)) %>%
  filter(enzyme_type == "kinase") %>%
  filter(translation %in% c("induced", "inhibited")) %>%
  select(GENE, translation) %>%
  table() %>%
  as.matrix()
kinase_order_tab <- kinase_order_tab[order(kinase_order_tab[,2], kinase_order_tab[,1], decreasing = T),]

tab2p <- kinase_translation_tab
tab2p <- tab2p %>%
  mutate(y = GENE) %>%
  mutate(x = translation) %>%
  mutate(point_size = Freq) %>%
  mutate(point_fill = translation) %>%
  filter(Freq > 0) %>%
  filter(x %in% c("induced", "inhibited"))
tab2p$y <- factor(tab2p$y, levels = rownames(kinase_order_tab))
tab2p_kinase <- tab2p

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

fn <- paste0(makeOutDir(resultD = resultD), "num_translation_pairs.pdf")
pdf(file = fn, width = 2.5, height = 10)
p
dev.off()

tab2p <- regression
tab2p <- tab2p %>%
  filter(regulated == T) %>%
  filter(!is.na(translation)) %>%
  filter(GENE %in% tab2p_kinase$GENE) %>%
  filter(translation %in% c("induced", "inhibited")) %>%
  mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
  select(GENE, phosphosite, translation, pair) %>%
  unique

tab2p <- melt(tab2p, id.vars = c("pair", "translation"))
tab2p <- tab2p %>%
  mutate(x = ifelse(grepl(pattern = "_", x = value), "phosphosite", "kinase")) %>%
  arrange(translation) %>%
  mutate(GENE = str_split_fixed(string = pair, pattern = ":", n = 3)[,1])
tab2p <- merge(tab2p, tab4lineplot_y, by = c("GENE"), all.x = T)


kinase_y <- 0
phosphosite_y <- 0
tab2p_new <- NULL
for (kinase_tmp in rownames(kinase_order_tab)) {
  tab2p_kinase_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "GENE")
  
  tab2p_phosphosite_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "phosphosite")
  
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

p <- ggplot()
p <- p + geom_text(data = tab2p[tab2p$x == "kinase",], mapping = aes(x = x, y = y, label = value), nudge_x = -0.2)
p <- p + geom_text(data = tab2p[tab2p$x == "phosphosite",], mapping = aes(x = x, y = y, label = value, color = translation), 
                   nudge_x = 0.3)
p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = pair, color = translation), alpha = 0.5)
p <- p + scale_color_manual(values = c("FALSE" =  "orange", "TRUE" = "grey50", "induced" = set1[1], "inhibited" = set1[2]))
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank(), 
               axis.title.y = element_blank(), 
               axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(), 
               axis.text.y = element_blank(), 
               axis.ticks.y = element_blank())
fn <- paste0(makeOutDir(resultD = resultD), "translation_pairs.pdf")
plot.new()
pdf(file = fn, width = 6, height = 10)
p
dev.off()













# get table and make figure about activity -----------------------------------------------
regression %>%
  select(ON_FUNCTION) %>%
  filter(grepl(pattern = "activity", x = ON_FUNCTION)) %>%
  unique

regression$activity <- sapply(1:nrow(regression), FUN = function(i, vector_function) {
  if (!is.na(vector_function[i]) & grepl(x = vector_function[i], pattern = "activity")) {
    tmp <- "altered"
    if (!grepl(x = vector_function[i], pattern = "activity, inhibited") & grepl(x = vector_function[i], pattern = "activity, induced")) {
      tmp <- "induced"
    }
    if (grepl(x = vector_function[i], pattern = "activity, inhibited") & !grepl(x = vector_function[i], pattern = "activity, induced")) {
      tmp <- "inhibited"
    }
  } else {
    tmp <- NA
  }
  return(tmp)
}, vector_function = regression$ON_FUNCTION)
kinase_activity_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(activity)) %>%
  filter(enzyme_type == "kinase") %>%
  select(GENE, activity) %>%
  table() %>%
  as.data.frame()
kinase_activity_tab


kinase_order_tab <- regression %>%
  filter(regulated == T) %>%
  filter(!is.na(activity)) %>%
  filter(enzyme_type == "kinase") %>%
  filter(activity %in% c("induced", "inhibited")) %>%
  select(GENE, activity) %>%
  table() %>%
  as.matrix() 
kinase_order_tab <- kinase_order_tab[order(kinase_order_tab[,2], kinase_order_tab[,1], decreasing = T),]

tab2p <- kinase_activity_tab
tab2p <- tab2p %>%
  mutate(y = GENE) %>%
  mutate(x = activity) %>%
  mutate(point_size = Freq) %>%
  mutate(point_fill = activity) %>%
  filter(Freq > 0) %>%
  filter(x %in% c("induced", "inhibited"))
tab2p$y <- factor(tab2p$y, levels = rownames(kinase_order_tab))
tab2p_kinase <- tab2p

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

fn <- paste0(makeOutDir(resultD = resultD), "num_activity_pairs.pdf")
pdf(file = fn, width = 2.5, height = 10)
p
dev.off()

tab2p <- regression
tab2p <- tab2p %>%
  filter(regulated == T) %>%
  filter(!is.na(activity)) %>%
  filter(GENE %in% tab2p_kinase$GENE) %>%
  filter(activity %in% c("induced", "inhibited")) %>%
  mutate(phosphosite = paste0(SUB_GENE, "_", SUB_MOD_RSD)) %>%
  select(GENE, phosphosite, activity, pair) %>%
  unique

tab2p <- melt(tab2p, id.vars = c("pair", "activity"))
tab2p <- tab2p %>%
  mutate(x = ifelse(grepl(pattern = "_", x = value), "phosphosite", "kinase")) %>%
  arrange(activity) %>%
  mutate(GENE = str_split_fixed(string = pair, pattern = ":", n = 3)[,1])
tab2p <- merge(tab2p, tab4lineplot_y, by = c("GENE"), all.x = T)


kinase_y <- 0
phosphosite_y <- 0
tab2p_new <- NULL
for (kinase_tmp in rownames(kinase_order_tab)) {
  tab2p_kinase_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "GENE")
  
  tab2p_phosphosite_tmp <- tab2p %>%
    filter(GENE == kinase_tmp) %>%
    filter(variable == "phosphosite")
  
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

p <- ggplot()
p <- p + geom_text(data = tab2p[tab2p$x == "kinase",], mapping = aes(x = x, y = y, label = value), nudge_x = -0.2)
p <- p + geom_text(data = tab2p[tab2p$x == "phosphosite",], mapping = aes(x = x, y = y, label = value, color = activity), 
                   nudge_x = 0.3)
p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = pair, color = activity, linetype = is.direct), alpha = 0.5)
p <- p + scale_color_manual(values = c("FALSE" =  "orange", "TRUE" = "grey50", "induced" = set1[1], "inhibited" = set1[2]))
p <- p + theme_bw()
p <- p + theme(axis.title.x = element_blank(), 
               axis.title.y = element_blank(), 
               axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(), 
               axis.text.y = element_blank(), 
               axis.ticks.y = element_blank())
fn <- paste0(makeOutDir(resultD = resultD), "activity_pairs.pdf")
pdf(file = fn, width = 6, height = 10)
p
dev.off()










