# Yige Wu @ WashU 2018 Jan
# make barplot showing kinase/phosphotase pairs with consistent phosphoprotein expression

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')

# input -----------------------------------------------------------------
us_diffexp <- fread(input = "cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/union_diffexp_regresson_genoalt/kinase_phosphotase_substrate_diffexp.txt")
us_diffexp$upstream_direction <- sapply(1:nrow(us_diffexp), function(n) ifelse(is.na(us_diffexp$FC_kin[n]), NA, ifelse(us_diffexp$FC_kin[n] > 1, "up", "down")))
us_diffexp$downstream_direction <- sapply(1:nrow(us_diffexp), function(n) ifelse(is.na(us_diffexp$FC_sub[n]), NA, ifelse(us_diffexp$FC_sub[n] > 1, "up", "down")))
us_diffexp$Cancer <- factor(us_diffexp$Cancer, levels = c("BRCA", "OV", "CO"))

# plot --------------------------------------------------------------------
us_diffexp_nona <- us_diffexp[!is.na(us_diffexp$upstream_direction) & !is.na(us_diffexp$downstream_direction),]
resultDnow <- makeOutDir()
dir.create(paste0(resultDnow, "free_scale"))
for (upstream in c("kinase", "phosphotase")) {
  for (upstream_direction in c("up", "down")) {
    for (downstream_direction in c("up", "down")) {
      df <- us_diffexp_nona[us_diffexp_nona$upstream == upstream & us_diffexp_nona$upstream_direction == upstream_direction & us_diffexp_nona$downstream_direction == downstream_direction, ]
      if (nrow(df) > 0) {
        figtitle <- paste(upstream, upstream_direction, "substrate", downstream_direction, "pairs", sep = "-")
        kinaseinpairs <- data.frame(table(df$KINASE))
        df <- merge(df, kinaseinpairs, by.x = c("KINASE"), by.y = c("Var1"), all.x = T)
        df$kinase <- reorder(df$KINASE, df$Freq)
        p <- ggplot(data = df, mapping = aes(x = kinase, fill = SUBSTRATE))
        p <- p + geom_bar()
        p <- p + coord_flip()
        p <- p + facet_wrap(~Cancer, scales = "free_x")
        p <- p + theme_nogrid()+ theme(legend.position="none") 
        p <- p + ggtitle(label = figtitle)
        p <- p + xlab(upstream) + ylab(paste0("count of substrate phosphosites"))
        p
        ggsave(filename = paste0(resultDnow, "free_scale/", figtitle, ".pdf"), height = 4, width = 6)
        
        p <- ggplot(data = df, mapping = aes(x = kinase, fill = SUBSTRATE))
        p <- p + geom_bar()
        p <- p + coord_flip()
        p <- p + facet_wrap(~Cancer)
        p <- p + theme_nogrid()+ theme(legend.position="none") 
        p <- p + ggtitle(label = figtitle)
        p <- p + xlab(upstream) + ylab(paste0("count of substrate phosphosites"))
        p
        ggsave(filename = paste0(resultDnow, figtitle, ".pdf"), height = 4, width = 6)
      }
    }
  }
}
