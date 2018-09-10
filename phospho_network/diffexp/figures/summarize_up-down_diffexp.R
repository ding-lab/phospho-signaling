# Yige Wu @ WashU 2018 Feb
# have a look at the how many proteins have both up/down-regulated phophosites

# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()
source(paste0(baseD, "cptac2p_analysis/pan3can_aes.R"))
Pho.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression/", "Pho.diffexp3can.txt"),
                         data.table = F)
library(dplyr)

# barplot showing the number/ratio of differentially phosphorylate --------
## initiate summary table
sum_table <- data.frame(table(Pho.diffexp3can[, c("diffexp_type", "Cancer")]))
for (cancer in c("BRCA", "OV", "CO")) {
  ## input protein and phosphoprotein data
  Pho <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_replicate_averaged_imputed.txt"),
               data.table = F)
  tmp <- data.frame(diffexp_type = "undiff", 
                    Cancer = cancer,
                    Freq = (nrow(Pho) - sum(sum_table$Freq[sum_table$Cancer == cancer])))
  sum_table <- rbind(sum_table, tmp)
}

sum_table$diffexp_type <- factor(sum_table$diffexp_type, levels = c("undiff", "up", "down"))
all_sites = sum_table %>%
  group_by(Cancer) %>%
  summarise(all_sites = sum(Freq))
all_sites <- data.frame(all_sites)
sum_table <- merge(sum_table, all_sites, all.x =T)
sum_table$ratio <- sum_table$Freq/sum_table$all_sites
sum_table$ratio_p <- paste0(sum_table$Freq ,"(", signif(100*sum_table$ratio, digits = 2), "%)")

p <- ggplot(data=sum_table, aes(y = Freq, x = Cancer, fill = diffexp_type, label = ratio_p))
p <- p + geom_bar(stat="identity", position='stack')
p <- p + geom_text(position = position_stack(vjust = 0.5))
p <- p + scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8", "undiff" = "grey"))
p <- p + theme_bw() + theme_nogrid() + scale_y_log10()
p <- p + ylab("number of differentially phosphorylated sites") +
  ggtitle("differentially phophorylated sites (FDR < 0.05, >2 fold change)")
p <- p + theme(axis.title.x = element_blank())
p <- p + theme(axis.text.x = element_text(colour="black", size=10, vjust=0.5), 
               axis.text.y = element_text(colour="black", size=10))
p
ggsave(filename = paste0(resultDnow, "cptac2p_3can_diffexp_phosphosites_up_down.pdf"),
       width = 6, height = 6)


# summarise down-regulated TSG --------------------------------------------


