# Yige Wu @ WashU Mar. 2019
## make pie charts showing the number of tumor samples


# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# Pie Chart with Percentages
slices <- c(123, 83, 97, 96, 105) 
# lbls <- c("breast cancer", "ovarian cancer", "colorectal cancer", "endometrial cancer", "clear cell\nrenal carcinoma")
lbls <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")
lbls <- paste0(lbls , "\n(n=", slices, ")") 
pdf(file = paste0(makeOutDir(resultD = resultD), "pie_number_samples.pdf"), width = 5, height = 5)
pie(slices, labels = lbls, col=rainbow(length(lbls)),
    main="number of tumor samples")
dev.off()