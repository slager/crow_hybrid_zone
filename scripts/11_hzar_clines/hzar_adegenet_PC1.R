## RUN THE TWO SCRIPTS IN THIS ORDER TO MAKE mtDNA nuDNA comparison GRAPH:
# 1) hzar_mtDNA.R
# 2) hzar_adegenet_PC1.R

library(adegenet)
library(hzar)
library(magrittr)

read.structure("Cb.s62.r0.05.maf0.02.het0.5.str.tsv.coastal.unlinked.nomd.str",
               n.ind=48,n.loc=905,onerowperind=F,
               col.lab=1,col.pop=0,col.others=0,
               row.marknames=1,NA.char="0") -> str

pop <- substr(rownames(str@tab),1,nchar(rownames(str@tab))-2)
#names(pop) <- NULL
str$pop <- as.factor(pop)
 
# ## PCA Coastal only
 
sum(is.na(str$tab))
X <- scaleGen(str,NA.method='mean')
pca1 <- dudi.pca(X,cent=F,scale=F,scannf=F,nf=2)
#barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
# # 
## Extract PCA1 scores
pca1$li$Axis1 -> p1
# flip so that 0=AK 1=CA
-p1 -> p1
p1
#Convert to (0,1) range standardized NOW FOR POPULATION MEANS
#Using linear transformation
#round(.0232921*p1+.415035,4) -> p

min0 =  -16.5497709  # mean of untransformed Homer -pc1
minf =  0

max0 =  20.4010785 #mean of untransformed ca -pc1
maxf =  1

a = (maxf-minf)/(max0-min0)
b = maxf-a*max0

p = a*p1+b

## Hybrid zone cline analysis

## Read distances for populations
read.csv("crow_geo_COASTAL_distance.csv",header=T,stringsAsFactors=F) -> d
## Population names
gsub("[0-9]","",rownames(str@tab)) -> popids
#unique(popids) -> popids
## Output vector of distances
d$distance -> dists
## Get population means
sapply(unique(popids),function(x){mean(p[which(popids %in% x)])},USE.NAMES=F) %>% round(4) -> means
#unname(means) -> means
#qqnorm(means)
## Get population variances
sapply(unique(popids),function(x){var(p[which(popids %in% x)])},USE.NAMES=F) -> vars
#unname(vars) -> vars
## Get population counts
sapply(unique(popids),function(x){length(p[which(popids %in% x)])},USE.NAMES=F) -> counts
#unname(counts) -> counts
data.frame(popids=unique(popids),dists,means,vars,counts) -> m

mydata <- hzar.doNormalData1DPops(distance=dists,siteID=unique(popids),muObs=means,varObs=vars,nEff=counts)
clineModel <- hzar.makeCline1DNormal(data=mydata,tails="none")
clineModel <- hzar.model.addCenterRange(clineModel, 1000,6000)
clineModel <- hzar.model.addMaxWidth(meta.model=clineModel,maxValue=7000)
#fitRequest <- hzar.first.fitRequest.old.ML(obsData=mydata,model=clineModel,verbose=T)  #Didn't work
fitRequest <- hzar.first.fitRequest.gC(gModel=clineModel,obsData=mydata,verbose=F)
#fitRequest$mcmcParam$chainLength <- 1e6
#fitRequest$mcmcParam$burnin <- 1e5
#fitRequest$mcmcParam$verbosity <- 0
#myfit <- hzar.doFit(fitRequest)
myfitlist_pca1 <- hzar.chain.doSeq(hzar.request=fitRequest,count=3,collapse=F)

'#CC3311' -> red
'#0077BB' -> blue
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c(blue,red))

#This adds a column of color values
# based on the y values
pca1_colors <- rbPal(100)[as.numeric(cut(means,breaks = 100))]
mtDNA_colors <- rbPal(100)[as.numeric(cut(pops$mtDNA_A,breaks = 100))]



pdf("cline_both2.pdf",5.5,5.5,useDingbats=F)
par(mar=c(4.1,4.1,4.1,2.6))
hzar.plot.cline(myfitlist_pca1[[3]],xlab="Pacific coastline (km)",ylab="Population value",xlim=c(0,7000),pch=NA)
hzar.plot.cline(myfitlist_mtDNA[[3]],add=T,lty=2,pch=NA)
# Get rid of line past points
rect( max(dists), .95, 7100, 1.05, col='white', lty=0)
rect( min(dists), -.05, -100, .05, col='white', lty=0)
points(dists,means,pch=3,col=pca1_colors) #plot pca1 data
points(pops$dist,pops$mtDNA_A,pch=1,col=mtDNA_colors) #plot mtDNA data
abline(v=c(3542,4864,5185,5656),lty=3,col='darkgray')
mtext(at=c(1671,4203,5024,5420,6428),text=c('AK','BC','A','OR','CA'),side=1,line=-1,cex=0.8)
mtext(at=c(5024),text=c('W'),side=1,line=-1.7,cex=0.8)
#mtext(at=c(0,7000),text=c('AK','CA'),side=1,line=-1)
legend(x=0,y=1.15,c("mtDNA fr(American)","Genomic mean PC1"),lty=c(2,1),pch=c(1,3),bg="white",cex=.9,bty="n")
dev.off()

# AK-BC  3542
# BC-WA  4864
# WA-OR  5185
# OR-CA  5656
# 
# 1771 AK
# 4203 BC
# 5024 WA
# 5420 OR
# 6328 CA
# =

#Fit to data group to enable next steps
hzar.fit2DataGroup(myfitlist_pca1[[3]]) -> fit3_pca1

#PDF graph version
pdf("pca1_cline_cred.pdf",6,6)
hzar.plot.fzCline(fit3_pca1,xlab="Smoothed coastline distance in km (AK to CA)",ylab="adegenet PC1",main="hzar.plot.cline()")
#hzar.plot.cline(myfitlist_pca1[[3]],xlab="Smoothed coastline distance in km (AK to CA)",ylab="adegenet PC1",main="hzar.plot.cline()")
dev.off()

hzar.get.ML.cline(myfitlist_pca1[[3]])$param.all$width
hzar.get.ML.cline(myfitlist_pca1[[3]])$param.all$center

#Extract center & +/- 2 LL range
hzar.getLLCutParam(fit3_pca1,"width",2)
hzar.getLLCutParam(fit3_pca1,"center",2)