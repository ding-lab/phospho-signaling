# source ------------------------------------------------------------------
packages = c(
  "ggplot2",
  "ggrepel",
  "ComplexHeatmap",
  "circlize",
  "RColorBrewer",
  "ggthemes",
  "rcartocolor",
  "Polychrome"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}


# color palettes ----------------------------------------------------------
RdBu = brewer.pal(9, "RdBu") 
RdBu1024 = colorRampPalette(rev(RdBu))(1024)


