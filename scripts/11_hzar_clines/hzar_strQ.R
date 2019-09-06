## RUN THE TWO SCRIPTS IN THIS ORDER TO MAKE mtDNA nuDNA comparison GRAPH:
# 1) hzar_mtDNA.R
# 2) hzar_strQ.R

library(adegenet)
library(hzar)
library(magrittr)
library(plyr)
library(dplyr)

## Get Structure K=2 values
df <- read.csv("parental define from Clumpp 0.98.csv",stringsAsFactors=F)
# Set population names
df$pop <- substr(df$id,1,nchar(df$id)-2)
# Choose Pacific coastal populations
df <- df %>%
  filter(pop %in% c("ca","cbc","ghc","hmr","jun","kit","nbc","neah","nvi","sea","vic","yvr"))
popids <- substr(df$id,1,nchar(df$id)-2)


df$A -> p
 
## Hybrid zone cline analysis

## Read distances for populations
read.csv("crow_geo_COASTAL_distance.csv",header=T,stringsAsFactors=F) -> d

## Output vector of distances
d$distance -> dists

## Get population means
sapply(unique(popids),function(x){mean(p[which(popids %in% x)])},USE.NAMES=F) %>% round(4) -> means

## Get population variances
sapply(unique(popids),function(x){var(p[which(popids %in% x)])},USE.NAMES=F) -> vars

## Get population counts
sapply(unique(popids),function(x){length(p[which(popids %in% x)])},USE.NAMES=F) -> counts

## Create data frame for input into Hzar
data.frame(popids=unique(popids),dists,means,vars,counts) -> m


mydata <- hzar.doNormalData1DPops(distance=dists,siteID=unique(popids),muObs=means,varObs=vars,nEff=counts)
clineModel <- hzar.makeCline1DNormal(data=mydata,tails="none")
clineModel <- hzar.model.addCenterRange(clineModel, 1000,6000)
clineModel <- hzar.model.addMaxWidth(meta.model=clineModel,maxValue=7000)

#Set initial mus based on observed terminus populations
mydata$frame['hmr','mu'] -> hzar.meta.init(clineModel)$muL
#0 -> hzar.meta.init(clineModel)$muL
mydata$frame['ca','mu'] -> hzar.meta.init(clineModel)$muR
#1 -> hzar.meta.init(clineModel)$muR

#Set initial variances based on observed terminus populations
mydata$frame['hmr','var'] -> hzar.meta.init(clineModel)$varL
mydata$frame['ca','var'] -> hzar.meta.init(clineModel)$varR

# Fix endpoint mus
hzar.meta.fix(clineModel)$muL <- TRUE
hzar.meta.fix(clineModel)$muR <- TRUE

fitRequest <- hzar.first.fitRequest.gC(gModel=clineModel,obsData=mydata,verbose=F)
myfitlist_strQ <- hzar.chain.doSeq(hzar.request=fitRequest,count=3,collapse=F)


'#CC3311' -> red
'#0077BB' -> blue
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c(blue,red))

#This adds a column of color values
# based on the y values
strQ_colors <- rbPal(100)[as.numeric(cut(means,breaks = 100))]
mtDNA_colors <- rbPal(100)[as.numeric(cut(pops$mtDNA_A,breaks = 100))]



pdf("cline_both2.pdf",5.5,5.5,useDingbats=F)
par(mar=c(4.1,4.1,4.1,2.6))
hzar.plot.cline(myfitlist_strQ[[3]],xlab="Pacific coastline (km)",ylab="Population value",xlim=c(0,7000),pch=NA,add=F)
hzar.plot.cline(myfitlist_mtDNA[[3]],add=T,lty=2,pch=NA)
# Get rid of line past points
rect( max(dists), .95, 7100, 1.05, col='white', lty=0)
rect( min(dists), -.05, -100, .05, col='white', lty=0)
points(dists,means,pch=3,col=strQ_colors) #plot strQ data
points(pops$dist,pops$mtDNA_A,pch=1,col=mtDNA_colors) #plot mtDNA data
abline(v=c(3542,4864,5185,5656),lty=3,col='darkgray')
mtext(at=c(1671,4203,5024,5420,6428),text=c('AK','BC','A','OR','CA'),side=1,line=-1,cex=0.8)
mtext(at=c(5024),text=c('W'),side=1,line=-1.7,cex=0.8)
#mtext(at=c(0,7000),text=c('AK','CA'),side=1,line=-1)
legend(x=0,y=1.1,c("mtDNA fr(American)","Genomic ancestry"),lty=c(2,1),pch=c(1,3),bg="white",cex=.9,bty="n")
dev.off()

#myfitlist_strQ[[3]]$mcmcRaw %>% tail

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
hzar.fit2DataGroup(myfitlist_strQ[[3]]) -> fit3_strQ

#PDF graph version
#pdf("strQ_cline_cred.pdf",6,6)
#hzar.plot.fzCline(fit3_strQ,xlab="Smoothed coastline distance in km (AK to CA)",ylab="adegenet PC1",main="hzar.plot.cline()")
#hzar.plot.cline(myfitlist_strQ[[3]],xlab="Smoothed coastline distance in km (AK to CA)",ylab="adegenet PC1",main="hzar.plot.cline()")
#dev.off()

hzar.get.ML.cline(myfitlist_strQ[[3]])$param.all$width
hzar.get.ML.cline(myfitlist_strQ[[3]])$param.all$center

#Extract center & +/- 2 LL range
hzar.getLLCutParam(fit3_strQ,"width",2)
hzar.getLLCutParam(fit3_strQ,"center",2)