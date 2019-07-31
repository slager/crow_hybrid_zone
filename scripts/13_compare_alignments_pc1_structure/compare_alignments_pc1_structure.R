library(adegenet)
library(magrittr)

## (1) Make data frame (df1) of transformed Genomic PC1 scores
## This is the coastal-only, no missing data alignment with 905 unlinked SNPs

read.structure("Cb.s62.r0.05.maf0.02.het0.5.str.tsv.coastal.unlinked.nomd.str",
               n.ind=48,n.loc=905,onerowperind=F,
               col.lab=1,col.pop=0,col.others=0,
               row.marknames=1,NA.char="0") -> str

pop <- substr(rownames(str@tab),1,nchar(rownames(str@tab))-2)
#names(pop) <- NULL
str$pop <- as.factor(pop)

# PCA Coastal only

sum(is.na(str$tab))
X <- scaleGen(str,NA.method='mean')
pca1 <- dudi.pca(X,cent=F,scale=F,scannf=F,nf=2)
#barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
# # 
## Extract PCA1 scores
pca1$li$Axis1 -> p1
# flip so that 0=AK 1=CA
-p1 -> p1

#Convert to (0,1) range standardized NOW FOR POPULATION MEANS
#Using linear transformation
#round(.0232921*p1+.415035,4) -> p

min0 =  -16.5497709  # mean of untransformed Homer -pc1
minf =  0

max0 =  20.4010785 #mean of untransformed ca -pc1
maxf =  1

a = (maxf-minf)/(max0-min0)
b = maxf-a*max0

p = a*p1+b  # p = transformed Genomic PC1 scores

# Get sample IDs in proper order from the adegenet object
attr(str@tab, "dimnames")[[1]] -> pca1_id # PC1 sample IDs

data.frame(id=pca1_id,adj_pca1=p) -> df1





## (2) Make data frame (df2) of Structure K=2 ancestry proportions.
## This is the alignment of 7,292 unlinked SNPs with coverage of at least 4/62 individuals continent-wide.
## This includes a much larger amount of missing data than the above alignment used for PCA

read.csv("parental define from Clumpp 0.98.csv") -> a

data.frame(id=a$id,structureA=a$A) -> df2

## (3) Merge df1, df2 for comparing Structure K=2 ancestry proportions and transformed Genomic PC1.
# Comapre the 48 Pacific Coastal samples
merge(df1,df2,by="id",all.x=TRUE) -> df

## (4) Plot Genomic PC1 vs. K=2 structure proportions
# and do regression

pdf("compare.pdf",6,6)
plot(df$adj_pca1,df$structureA,
     xlab="Genomic PC1",
     ylab="Structure K=2 ancestry")
abline(a=0,b=1)
dev.off()
summary(lm(df$structureA ~ df$adj_pca1))

# For regression, Adj. R squared is 0.9478, p < 2.2e-16

## Adjusted PC1 is from the 905-SNP alignment with no missing data
## Structure K=2 American ancestry proportion is from 7,292-SNP alignment with missing data for up to 58/62 contintent-wide samples
## This is for the 48 coastal samples