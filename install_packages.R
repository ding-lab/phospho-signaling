install.packages("stringr")
install.packages("reshape", dependencies=TRUE)
install.packages("data.table", dependencies=TRUE)
install.packages("readr", dependencies = T)
install.packages("gplots", dependencies = T)
install.packages("GO.db")
install.packages("eulerr")
install.packages("readr")
install.packages("readxl")
install.packages("impute")
install.packages("WGCNA")
install.packages("ggrepel")

source("https://bioconductor.org/biocLite.R")
biocLite("GO.db")
biocLite("impute")
biocLite("preprocessCore")

for (pkg_name_tmp in c("survival", "survminer", "cmprsk", "data.table", "pheatmap")) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp)
  }
}



