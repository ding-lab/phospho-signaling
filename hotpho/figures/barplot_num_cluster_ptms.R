# Yige Wu @ WashU 2018 Aug
## plot per cancer # hotpho clusters and number of ptms

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")

# inputs hotpho results (hg38 based)------------------------------------------------------------------
hotpho <- read_delim("~/Box Sync/Ding_Lab/Projects_Current/hotpho_data/output/Data_201807_cc.p0.05.cluster_transcriptSynced_hybrid_SourceAnnotated.tsv", 
                     "\t", escape_double = FALSE, col_types = cols(Cluster = col_character()), trim_ws = TRUE)

driver_table <- read_excel("./Ding_Lab/Projects_Current/TCGA_data/gene_lists/mmc1.xlsx", 
                           sheet = "Table S1", skip = 3)
driver_table <- data.frame(driver_table)


# take hotpho clusters with ptms annotated to CPTAC dataset and split per cohort---------------
hotpho_cptac <- hotpho[!is.na(hotpho$Cohort) & hotpho$Cohort != "Not_Found",]
table(hotpho_cptac$Cohort)
cohort_list <- sapply(hotpho_cptac$Cohort, FUN = function(c) strsplit(x = c, split = "\\|")[[1]], simplify = T)
cohort <- unlist(cohort_list, use.names = F)
cohort_num <- sapply(cohort_list, FUN = function(c) length(c))
hotpho_cptac_m <- hotpho_cptac[rep(x = 1:nrow(hotpho_cptac), cohort_num),]
hotpho_cptac_m$Cohort <- cohort
## remove the retrospective ones
hotpho_cptac_m <- hotpho_cptac_m[!grepl(pattern = "Retrospective", x = hotpho_cptac_m$Cohort) & !grepl(pattern = "CPTAC3", x = hotpho_cptac_m$Cohort),]
hotpho_cptac_m$Cohort <- as.vector(str_split_fixed(string = hotpho_cptac_m$Cohort, pattern = "-", n = 2)[, 2])
hotpho_cptac_m$rsd <- str_split_fixed(string = hotpho_cptac_m$Mutation_Gene, pattern = "\\.", 2)[, 2]

length(unique(hotpho_cptac_m$Gene_Drug))
# set variables -----------------------------------------------------------
## assign colors to each cohorts
color_cohorts <- brewer.pal(length(unique(hotpho_cptac_m$Cohort)), "Set1")
names(color_cohorts) <- unique(hotpho_cptac_m$Cohort)
top_genes2show <- 20

# plot # clusters (y) per cohort (x)----------------------------------------
count_tab <- data.frame(table(unique(hotpho_cptac_m[, c("Cohort", "Cluster")])[, "Cohort"]))
tab2p <- count_tab
tab2p$Var1 <- factor(tab2p$Var1, levels = as.vector(tab2p$Var1)[order(tab2p$Freq, decreasing = T)])
p <- ggplot()
p <- p + geom_bar(data=tab2p, aes(y = Freq, x = Var1, fill = Var1), 
                  stat="identity", position='stack', color = NA)
p <- p + scale_fill_manual(values = color_cohorts)
p <- p + geom_text(data = tab2p, mapping = aes(x = Var1, y = Freq, label = Freq), nudge_y = 4, size = 5, color = "black")
p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
p <- p + xlab("Cohort")+ylab("No.hybrid_clusters")
p <- p + theme(axis.title=element_text(size=15, face = "bold"))
p <- p + theme(axis.text.x = element_blank())
p <- p + theme(axis.text.y = element_text(colour="black", size=10))
p
fn = paste0(makeOutDir(resultD = resultD), 'No.hybrid_clusters_per_Cohort.pdf')
ggsave(file=fn, height=4, width=7)

# plot # clusters (y) per top gene (x) per cancer type (color)----------------------------------------
count_tab4order <- data.frame(table(unique(hotpho_cptac_m[, c("Cohort", "Cluster", "Gene_Drug")])[, c("Gene_Drug")]))
count_tab4order <- count_tab4order[order(count_tab4order$Freq, decreasing = T), ]
top_x <- (count_tab4order$Var1[1:top_genes2show])

count_tab <- data.frame(table(unique(hotpho_cptac_m[, c("Cohort", "Cluster", "Gene_Drug")])[, c("Cohort", "Gene_Drug")]))
count_tab <- count_tab[count_tab$Gene_Drug %in% top_x,]

tab2p <- count_tab
tab2p$Gene_Drug <- factor(tab2p$Gene_Drug, levels = top_x)
p <- ggplot()
p <- p + geom_bar(data=tab2p, aes(y = Freq, x = Gene_Drug, fill = Cohort), 
                  stat="identity", position='stack', color = NA)
p <- p + scale_fill_manual(values = color_cohorts)
# p <- p + geom_text(data = tab2p, mapping = aes(x = Var1, y = Freq, label = Freq), nudge_y = 4, size = 5, color = "black")
p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
p <- p + xlab("Cohort")+ylab("No.hybrid_clusters")
p <- p + theme(axis.title.y=element_text(size=15, face = "bold"))
p <- p + theme(axis.title.x=element_blank())
p <- p + theme(axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 0.5, angle = 90, face = "bold"))
p <- p + theme(axis.text.y = element_text(colour="black", size=10))
p
fn = paste0(makeOutDir(resultD = resultD), 'No.hybrid_clusters_per_top',top_genes2show,'_hybrid_clusters_gene_per_Cohort.pdf')
ggsave(file=fn, height=4, width=7)

count_tab <- data.frame(table(unique(hotpho_cptac_m[, c("Cohort", "Position", "Gene_Drug")])[, c("Cohort", "Gene_Drug")]))
count_tab <- count_tab[count_tab$Gene_Drug %in% top_x,]

tab2p <- count_tab
tab2p$Gene_Drug <- factor(tab2p$Gene_Drug, levels = top_x)
p <- ggplot()
p <- p + geom_bar(data=tab2p, aes(y = Freq, x = Gene_Drug, fill = Cohort), 
                  stat="identity", position='stack', color = NA)
p <- p + scale_fill_manual(values = color_cohorts)
# p <- p + geom_text(data = tab2p, mapping = aes(x = Var1, y = Freq, label = Freq), nudge_y = 4, size = 5, color = "black")
p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
p <- p + xlab("Cohort")+ylab("No.PTM_sites_in_hybrid_clusters")
p <- p + theme(axis.title.y=element_text(size=10, face = "bold"))
p <- p + theme(axis.title.x=element_blank())
p <- p + theme(axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 0.5, angle = 90, face = "bold"))
p <- p + theme(axis.text.y = element_text(colour="black", size=10))
p
fn = paste0(makeOutDir(resultD = resultD), 'No.PTM_sites_in_hybrid_clusters_per_top',top_genes2show,'_hybrid_clusterses_gene_per_Cohort.pdf')
ggsave(file=fn, height=4, width=7)

# plot # ptms sites (y) per cohort (x) ----------------------------------------
count_tab <- data.frame(table(unique(hotpho_cptac_m[, c("Cohort", "Position")])[, "Cohort"]))
tab2p <- count_tab
tab2p$Var1 <- factor(tab2p$Var1, levels = as.vector(tab2p$Var1)[order(tab2p$Freq, decreasing = T)])
p <- ggplot()
p <- p + geom_bar(data=tab2p, aes(y = Freq, x = Var1, fill = Var1), 
                  stat="identity", position='stack', color = NA)
p <- p + scale_fill_manual(values = color_cohorts)
p <- p + geom_text(data = tab2p, mapping = aes(x = Var1, y = Freq, label = Freq), nudge_y = 4, size = 5, color = "black")
p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
p <- p + xlab("Cohort")+ylab("No.PTM_sites_in_hybrid_clusters")
p <- p + theme(axis.title=element_text(size=15, face = "bold"))
p <- p + theme(axis.text.x = element_blank())
p <- p + theme(axis.text.y = element_text(colour="black", size=10))
p
fn = paste0(makeOutDir(resultD = resultD), 'No.PTM_sites_in_hybrid_clusters_per_Cohort.pdf')
ggsave(file=fn, height=4, width=7)

# plot # ptms sites (y) per top gene (x) per cohort (color) ----------------------------------------
count_tab4order <- data.frame(table(unique(hotpho_cptac_m[, c("Cohort", "Position", "Gene_Drug")])[, c("Gene_Drug")]))
count_tab4order <- count_tab4order[order(count_tab4order$Freq, decreasing = T), ]
top_x <- (count_tab4order$Var1[1:top_genes2show])

count_tab <- data.frame(table(unique(hotpho_cptac_m[, c("Cohort", "Position", "Gene_Drug")])[, c("Cohort", "Gene_Drug")]))
count_tab <- count_tab[count_tab$Gene_Drug %in% top_x,]

tab2p <- count_tab
tab2p$Gene_Drug <- factor(tab2p$Gene_Drug, levels = top_x)
p <- ggplot()
p <- p + geom_bar(data=tab2p, aes(y = Freq, x = Gene_Drug, fill = Cohort), 
                  stat="identity", position='stack', color = NA)
p <- p + scale_fill_manual(values = color_cohorts)
# p <- p + geom_text(data = tab2p, mapping = aes(x = Var1, y = Freq, label = Freq), nudge_y = 4, size = 5, color = "black")
p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
p <- p + xlab("Cohort")+ylab("No.PTM_sites_in_hybrid_clusters")
p <- p + theme(axis.title.y=element_text(size=10, face = "bold"))
p <- p + theme(axis.title.x=element_blank())
p <- p + theme(axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 0.5, angle = 90, face = "bold"))
p <- p + theme(axis.text.y = element_text(colour="black", size=10))
p
fn = paste0(makeOutDir(resultD = resultD), 'No.PTM_sites_in_hybrid_clusters_per_top',top_genes2show,'_PTM_sites_gene_per_Cohort.pdf')
ggsave(file=fn, height=4, width=7)


count_tab <- data.frame(table(unique(hotpho_cptac_m[, c("Cohort", "Cluster", "Gene_Drug")])[, c("Cohort", "Gene_Drug")]))
count_tab <- count_tab[count_tab$Gene_Drug %in% top_x,]

tab2p <- count_tab
tab2p$Gene_Drug <- factor(tab2p$Gene_Drug, levels = top_x)
p <- ggplot()
p <- p + geom_bar(data=tab2p, aes(y = Freq, x = Gene_Drug, fill = Cohort), 
                  stat="identity", position='stack', color = NA)
p <- p + scale_fill_manual(values = color_cohorts)
# p <- p + geom_text(data = tab2p, mapping = aes(x = Var1, y = Freq, label = Freq), nudge_y = 4, size = 5, color = "black")
p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
p <- p + xlab("Cohort")+ylab("No.hybrid_clusters")
p <- p + theme(axis.title.y=element_text(size=15, face = "bold"))
p <- p + theme(axis.title.x=element_blank())
p <- p + theme(axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 0.5, angle = 90, face = "bold"))
p <- p + theme(axis.text.y = element_text(colour="black", size=10))
p
fn = paste0(makeOutDir(resultD = resultD), 'No.hybrid_clusters_per_top',top_genes2show,'_PTM_sites_gene_per_Cohort.pdf')
ggsave(file=fn, height=4, width=7)


# Bubble plots ------------------------------------------------------------
count_tab4orderx <- data.frame(table(unique(hotpho_cptac_m[hotpho_cptac_m$Gene_Drug %in% driver_table$Gene, c("Cohort", "Position", "Gene_Drug")])[, c("Gene_Drug")]))
count_tab4orderx <- count_tab4orderx[order(count_tab4orderx$Freq, decreasing = T), ]
top_x <- as.vector(count_tab4orderx$Var1)
count_tab4ordery <- data.frame(table(unique(hotpho_cptac_m[hotpho_cptac_m$Gene_Drug %in% driver_table$Gene, c("Cohort", "Position", "Gene_Drug")])[, c("Cohort")]))
count_tab4ordery <- count_tab4ordery[order(count_tab4ordery$Freq, decreasing = T), ]
top_y <- as.vector(count_tab4ordery$Var1)

count_tab <- data.frame(table(unique(hotpho_cptac_m[, c("Cohort", "Position", "Gene_Drug")])[, c("Cohort", "Gene_Drug")]))
count_tab <- count_tab[count_tab$Gene_Drug %in% top_x,]
count_tab <- count_tab[count_tab$Freq > 0,]
tab2p <- count_tab
tab2p$Gene_Drug <- factor(tab2p$Gene_Drug, levels = top_x)
tab2p$Cohort <- factor(tab2p$Cohort, levels = rev(top_y))

p <- ggplot()
p <- p + geom_point(data=tab2p, aes(y = Cohort, x = Gene_Drug, color = Cohort, size = Freq), alpha = 0.8)
p <- p + scale_color_manual(values = color_cohorts)
p <- p + geom_text(data = tab2p, mapping = aes(y = Cohort, x = Gene_Drug, label = Freq), size = 2, color = "black")
p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
p <- p + theme(axis.title.x=element_blank())
p <- p + theme(axis.title.y=element_blank())
p <- p + theme(axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 0.5, angle = 90, face = "bold"))
p <- p + theme(axis.text.y = element_text(colour="black", size=15, face = "bold"))
p
fn = paste0(makeOutDir(resultD = resultD), 'No.PTM_sites_in_hybrid_clusters_per_top',top_genes2show,'_PTM_sites_gene_per_Cohort_bubble.pdf')
ggsave(file=fn, height=3, width=4)



