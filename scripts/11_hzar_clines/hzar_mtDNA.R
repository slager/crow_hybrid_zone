## RUN THE TWO SCRIPTS IN THIS ORDER TO MAKE mtDNA nuDNA comparison GRAPH:
# 1) hzar_mtDNA.R
# 2) hzar_adegenet_PC1.R

library(hzar)

read.csv("pop_haps_pacific.csv",stringsAsFactors=F) -> pops

mydata <- hzar.doMolecularData1DPops(
  distance=pops$dist,
  pObs=pops$mtDNA_A,
  nEff=pops$n,
  siteID=pops$pop)

clineModel <- hzar.makeCline1DFreq(data=mydata,tails="none")
clineModel <- hzar.model.addCenterRange(clineModel, 1000,6000)
clineModel <- hzar.model.addMaxWidth(meta.model=clineModel,maxValue=7000)
fitRequest <- hzar.first.fitRequest.old.ML(model=clineModel,obsData=mydata,verbose=T)  #Didn't work
#fitRequest <- hzar.first.fitRequest.gC(gModel=clineModel,obsData=mydata,verbose=F)
#fitRequest$mcmcParam$chainLength <- 1e6
#fitRequest$mcmcParam$burnin <- 1e5
#fitRequest$mcmcParam$verbosity <- 0
#myfit <- hzar.doFit(fitRequest)
myfitlist_mtDNA <- hzar.chain.doSeq(hzar.request=fitRequest,count=3,collapse=F)
hzar.plot.cline(myfitlist_mtDNA[[3]],xlab="Smoothed coastline distance in km (AK to CA)",ylab="American mtDNA",main="hzar.plot.cline()")


#Fit to data group to enable next steps
hzar.fit2DataGroup(myfitlist_mtDNA[[3]]) -> fit3_mtDNA

#PDF graph version
pdf("cline_cred_mtDNA.pdf",6,6)
hzar.plot.fzCline(fit3_mtDNA,xlab="Smoothed coastline distance in km (AK to CA)",ylab="American mtDNA",main="hzar.plot.cline()")
#hzar.plot.cline(myfitlist_mtDNA[[3]],xlab="Smoothed coastline distance in km (AK to CA)",ylab="American mtDNA",main="hzar.plot.cline()")
dev.off()


hzar.get.ML.cline(myfitlist_mtDNA[[3]])$param.all$width
hzar.get.ML.cline(myfitlist_mtDNA[[3]])$param.all$center

#Extract center & +/- 2 LL range
hzar.getLLCutParam(fit3_mtDNA,"width",2)
hzar.getLLCutParam(fit3_mtDNA,"center",2)
