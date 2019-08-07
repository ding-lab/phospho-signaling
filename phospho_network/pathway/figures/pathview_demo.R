library(pathview)
data(gse16873.d)
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",
                   species = "hsa", out.suffix = "gse16873")
(pv.out)
