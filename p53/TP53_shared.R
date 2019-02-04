



# TP53 domain -------------------------------------------------------------


getTP53domain <- function(position) {
  if (is.na(position)) {
    return("other")
  } else if (position <= 43) {
    return("TAD1")
  } else if (position <= 63) {
    return("TAD2")
  } else if (position <= 92) {
    return("proline rich\ndomain")
  } else  if (position >= 102 & position <= 292) {
    return("DBD")
  } else if (position >= 316 & position <= 325) {
    return("NLS")
  } else if (position >= 325 & position <= 355) {
    return("tetramerization\ndomain")
  } else if (position >= 356 & position <= 393) {
    return("regulatory\ndomain")
  } else {
    return("other")
  }
}

getTP53domains <- function(positions) {
  domains <- sapply(positions, getTP53domain)
  return(domains)
}

# plotting ----------------------------------------------------------------
library(RColorBrewer)


if(!("pheatmap" %in% installed.packages()[,"Package"])) {
  install.packages("pheatmap")
}
library(pheatmap)
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf <- function(x, filename, width=6, height=6) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}