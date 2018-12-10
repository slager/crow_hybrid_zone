library(adegenet)


read.structure("Cb.s62.r0.05.maf0.02.het0.5.tsv.unlinked.nomd.str",
               n.ind=62,n.loc=738,onerowperind=F,
               col.lab=1,col.pop=0,col.others=0,
               row.marknames=1,NA.char="0") -> str

str(str)
class(str@tab)
str@tab[1:4,1:4]
rownames(str@tab)

grp <- find.clusters(str, max.n.clust=10,n.pca=1e6)  # Manually enter 2, press <Enter>

grp$Kstat

