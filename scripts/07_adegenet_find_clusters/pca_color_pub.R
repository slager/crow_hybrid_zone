library(adegenet)
#library(rgdal)
#library(sp)
#library(MASS)
#library(RColorBrewer)

read.structure("Cb.s62.r0.05.maf0.02.het0.5.tsv.unlinked.nomd.str",
               n.ind=62,n.loc=738,onerowperind=F,
               col.lab=1,col.pop=0,col.others=0,
               row.marknames=1,NA.char="0") -> str

str(str)
class(str@tab)
str@tab[1:4,1:4]
rownames(str@tab)

sum(is.na(str$tab))
X <- scaleGen(str,NA.method='mean')
pca1 <- dudi.pca(X,cent=F,scale=F,scannf=F,nf=2)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
# 
# pdf("PCA.labels.pdf",15,15)
# s.label(pca1$li,grid=F,boxes=T,addaxes=T)
# #title("PCA of microbov dataset\naxes 1-2")
# title("738 full-coverage SNPs",xlab="PCA1",ylab="PCA2")
# add.scatter.eig(pca1$eig[1:20], 2,1,2)
# dev.off()
# 
# pca1$eig
# pca1$li

pop <- substr(rownames(str@tab),1,nchar(rownames(str@tab))-2)
#names(pop) <- NULL
str$pop <- as.factor(pop)
# 
# ## Black-and-white
# pdf("PCA.pops.b-w.pdf",10,10)
# s.class(pca1$li, str$pop)  ## See IBD section below for str$pop
# #title("PCA of microbov dataset\naxes 1-2")
# title("738 full-coverage SNPs",xlab="PCA1",ylab="PCA2")
# add.scatter.eig(pca1$eig[1:20], 3,1,2)
# dev.off()

## Color
# 738 full coverage SNPs
pdf("PCA.pops.color.pdf",10,10)
#col <- heat.colors(length(levels(str$pop)))
#col <- rep(c('black','darkgreen','darkblue','darkmagenta','darkorange','firebrick','gold3','deeppink'),10)[sample(1:length(levels(str$pop)),length(levels(str$pop)),replace=F)]
col <- c("darkblue", #ca
         "deeppink", #cbc
         "gold3", #ewa
         "black", #ghc
         "darkorange", #hmr
         "firebrick", #jun
         "darkorange", #kit
         "black", #la
         "gold3", #mi
         "darkgreen", #nbc
         "darkblue",#neah
         "deeppink", #nf
         "black",#nvi
         "firebrick",#nynj
         "darkmagenta",#sea
         "darkgreen",#vic
         "gold3")#yvr
s.class(pca1$li,str$pop,col=transp(col,1),grid=F,addaxes=F)  ## See also IBD section for str$pop
#title("PCA of microbov dataset\naxes 1-2")
#title(xlab="PC1 (6.8%)",ylab="PC2 (3.3%)")
mtext("Genomic PC1 (6.8%)",side=1,line=3,cex=1.5)
mtext("Genomic PC2 (3.3%)",side=2,line=2,cex=1.5)
add.scatter.eig(w=pca1$eig,nf=2,xax=1,yax=2)
dev.off()

(pca1$eig/sum(pca1$eig))[1:2] # Variation explained by PC1,PC2
