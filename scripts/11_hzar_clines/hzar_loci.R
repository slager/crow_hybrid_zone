library(adegenet)
library(hzar)
library(parallel)
library(magrittr)
library(plyr)
library(dplyr)
library(scales)

read.structure("Cb.s62.r0.05.maf0.02.het0.5.str.tsv.coastal.unlinked.nomd.str",
               n.ind=48,n.loc=905,onerowperind=F,
               col.lab=1,col.pop=0,col.others=0,
               row.marknames=1,NA.char="0") -> str

pop <- substr(rownames(str@tab),1,nchar(rownames(str@tab))-2)
#names(pop) <- NULL
str$pop <- as.factor(pop)


## Hybrid zone cline analysis

## Read distances for populations
read.csv("crow_geo_COASTAL_distance.csv",header=T,stringsAsFactors=F) -> d
## Population names
gsub("[0-9]","",rownames(str@tab)) -> popids
#unique(popids) -> popids
## Output vector of distances
d$distance -> dists

####
####
####


## Organize empty list structures for loci

locus.names <- names(str@all.names)  ## Full dataset
#locus.names <- names(str@all.names)[1:2] ## For testing
loci <- vector("list", length(locus.names))
names(loci) <- locus.names
for (i in names(loci)){
  list() -> loci[[i]]$ref_allele
  list() -> loci[[i]]$alt_allele
  list() -> loci[[i]]$obs
  list() -> loci[[i]]$models
  list() -> loci[[i]]$fitRs
  list() -> loci[[i]]$runs
  list() -> loci[[i]]$analysis
}

## Create hzar.obsData object in loci[[i]]$obs
for (i in names(loci)){
  grablocus <- i
  grabcolnames <- colnames(str@tab)[which(gsub("[.][0-9]$","",colnames(str@tab))==grablocus)]
  graballele <- grabcolnames[1]  ## 1st allele is Reference Allele.
  loci[[i]]$ref_allele <- graballele
  loci[[i]]$alt_allele <- grabcolnames[2] ## 2nd allele is Alternate Allele
  for (p in unique(popids)){
    grabpop <- p
    grabrownames <- rownames(str@tab)[which(gsub("[0-9]","",rownames(str@tab))==grabpop)]
    grabbed <- str$tab[grabrownames,graballele]
    sum(grabbed)/(2*length(grabbed)) -> alleleFreq
    2*length(grabbed) -> alleleN
  }
  sapply(unique(popids),function(grabpop){
    grabrownames <- rownames(str@tab)[which(gsub("[0-9]","",rownames(str@tab))==grabpop)]
    grabbed <- str$tab[grabrownames,graballele]
    sum(grabbed)/(2*length(grabbed)) -> alleleFreq
    2*length(grabbed) -> alleleN
    return(c(alleleFreqs=alleleFreq,alleleNs=alleleN))
  }) -> f.n
  t(f.n)[,'alleleFreqs'] -> alleleFreqs
  t(f.n)[,'alleleNs'] -> alleleNs
  #print(length(alleleNs))
  #print(length(alleleFreqs))
  loci[[i]]$obs <- hzar.doMolecularData1DPops(distance=dists,pObs=alleleFreqs,nEff=alleleNs,siteID=unique(popids))
}

## Look at data
#for (i in 1:10){   ## Check out the 1st 50 loci
#hzar.plot.obsData(loci[[i]]$obs,main=loci[[i]]$ref_allele)
#}

##########
## TESTING:  Use subset of specific loci ##
##########

#loci[c("991_21","73_37","25_12")] -> loci

## Load the models to be run

# Helper function
loadClineModel <- function(locus,scaling,tails,id=paste(scaling,tails,sep=".")){
  loci[[locus]]$models[[id]] <<- hzar.makeCline1DFreq(data=loci[[locus]]$obs, scaling, tails,direction=NULL)
}

for (i in names(loci)){
  #loadClineModel(locus=i,scaling='fixed',tails='none')
  loadClineModel(locus=i,scaling='free',tails='none')
  #loadClineModel(locus=i,scaling='fixed',tails='mirror')
  #loadClineModel(locus=i,scaling='fixed',tails='left')
  #loadClineModel(locus=i,scaling='fixed',tails='right')
  #loadClineModel(locus=i,scaling='fixed',tails='both')
}

## Set boundaries for center and width parameters
for (i in names(loci)){
  loci[[i]]$models <- sapply(loci[[i]]$models,hzar.model.addBoxReq,0,4000,simplify=FALSE)
}

Sys.time() ## Starting cluster processes

## Compile each of the models to prepare for fitting. Takes longer, do in parallel
#  i.e., first fit request
#  Do in 2 steps to avoid scoping problems
#  Let parLapply do "reduce" step, then assign to main list in GlobalEnv afterwards

makeCluster(max(1,detectCores()-1),type="FORK") -> cluster ## Keeps active packages & variables set. More memory-efficient.
clusterExport(cluster,"loci")
globalList <- parLapply(cluster,loci,function(locus){
  sapply(locus$models,hzar.first.fitRequest.old.ML,obsData=locus$obs,verbose=FALSE,simplify=FALSE)
})
stopCluster(cluster)
for (i in names(loci)) {globalList[[i]] -> loci[[i]]$fitRs$init}


## Do initial runs of each model
for (i in names(loci)) {loci[[i]]$runs$init <- list() }
makeCluster(max(1,detectCores()-1),type="FORK") -> cluster ## Keeps active packages & variables set. More memory-efficient.
clusterExport(cluster,"loci")
globalList <- parLapply(cluster,loci,function(locus){
  sapply(locus$fitRs$init,hzar.doFit,simplify=FALSE)
})
stopCluster(cluster)
for (i in names(loci)) {globalList[[i]] -> loci[[i]]$runs$init}
rm(globalList)

#loci -> loci_backup
#myfitlist <- hzar.chain.doSeq(hzar.request=fitRequest,count=3,collapse=F)

### String 2 more dependent chains for each model, starting with initial run
for (i in names(loci)) {loci[[i]]$runs$chain <- list() }
makeCluster(max(1,detectCores()-1),type="FORK") -> cluster ## Keeps active packages & variables set. More memory-efficient.
clusterExport(cluster,"loci")
globalList <- parLapply(cluster,loci,function(locus){
  sapply(locus$runs$init,hzar.chain.doSeq,count=2,collapse=F,simplify=FALSE)
})
stopCluster(cluster)
for (i in names(loci)) {globalList[[i]] -> loci[[i]]$runs$chain}
rm(globalList)

Sys.time()  ## Stopping cluster processes

## Plot traces for 1st 3 loci
#for (i in 1:3){
#  plot(hzar.mcmc.bindLL(loci[[i]]$runs$chain$none[[2]]))
#}

## Plot clines for 1st 3 loci
#for (i in 1:3){
#  hzar.plot.cline(loci[[i]]$runs$chain$none[[2]],main=loci[[i]]$ref_allele)
#}

## Start aggregating the analyses for AIC

# Create a model data group for the null model (expected allele
# frequency independent of distance along cline) to include in
# analysis.
for (i in names(loci)) {loci[[i]]$analysis$initDGs <- list(nullModel = hzar.dataGroup.null(loci[[i]]$obs)) }

# Create a model data group (hzar.dataGroup object) for each
# model from the initial runs.

for (i in names(loci)) {
  for (model in names(loci[[i]]$runs$init)) {
    loci[[i]]$analysis$initDGs[[model]] <-hzar.dataGroup.add(loci[[i]]$runs$init[[model]])
  }
}

#loci -> loci_backup

## Create a hzar.obsDataGroup object from the four hzar.dataGroup
## just created, copying the naming scheme 
## LONG, consider making parallel
for (i in names(loci)) {loci[[i]]$analysis$oDG <-hzar.make.obsDataGroup(loci[[i]]$analysis$initDGs)}
for (i in names(loci)) {loci[[i]]$analysis$oDG <-  hzar.copyModelLabels(loci[[i]]$analysis$initDGs,loci[[i]]$analysis$oDG)}

## Convert all runs to hzar.dataGroup objects, adding them to
## the hzar.obsDataGroup object.
## LONG, consider making parallel
for (i in names(loci)){loci[[i]]$analysis$oDG <- hzar.make.obsDataGroup(lapply(loci[[i]]$runs$chain,hzar.dataGroup.add),loci[[i]]$analysis$oDG)}

## Check names of hzar.dataGroup objects in the hzar.obsDataGroup object.
print(summary(loci[[1]]$analysis$oDG$data.groups))

## Do model selection of cline/null models based on the AICc scores
for (i in names(loci)){
print(loci[[i]]$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(loci[[i]]$analysis$oDG))
}

## Use/Print out the model with the minimum AICc score
#for (i in names(loci)){
#print(loci[[i]]$analysis$model.name <- 
#        rownames(loci[[i]]$analysis$AICcTable)[[ which.min(loci[[i]]$analysis$AICcTable$AICc )]])
#}

## Conservatively use nullModel unless another model has dAICc > 6
for (i in names(loci)){
  loci[[i]]$analysis$AICcTable -> d
  d[order(d$AICc),,drop=FALSE] -> d
  d$AICc - min(d$AICc) -> d$dAICc
  if ("nullModel" %in% rownames(d[which(d$dAICc < 6),])) {"nullModel" -> m} else { rownames(d[which.min(d$dAICc),]) -> m }
  print(loci[[i]]$analysis$model.name <- m)
}

## Extract the hzar.dataGroup object for the selected model
for (i in names(loci)){
loci[[i]]$analysis$model.selected <- loci[[i]]$analysis$oDG$data.groups[[loci[[i]]$analysis$model.name]]
}

## Look at the variation in parameters for the selected model
#for (i in names(loci)){
#  if (loci[[i]]$analysis$model.name != "nullModel"){
#    print(hzar.getLLCutParam(loci[[i]]$analysis$model.selected,names(loci[[i]]$analysis$model.selected$data.param)))
#  }
#}

## Print the maximum likelihood cline for the selected model
#for (i in names(loci)){
#print(hzar.get.ML.cline(loci[[i]]$analysis$model.selected))
#}

## Get single Maximum Likelihood Parameters
#hzar.get.ML.cline(loci[[2]]$analysis$model.selected)$param.free
#hzar.get.ML.cline(loci[[4]]$analysis$model.selected)$param.all

Sys.time()  ## Before plotting




## Extract data frame of parameters
cline.data <-
data.frame(ref_allele=sapply(loci,function(i){i[['ref_allele']]
           }),
           alt_allele=sapply(loci,function(i){i[['alt_allele']]
           }),
           model=sapply(loci,function(i){
             i$analysis$model.name
           }),
           K=sapply(loci,function(i){
             length(hzar.get.ML.cline(i$analysis$model.selected)$param.all)
           }),
           center2LLLow=sapply(loci,function(i){
             ifelse(
               is.null(hzar.get.ML.cline(i$analysis$model.selected)$param.all$center),
               NA,
               hzar.getLLCutParam(i$analysis$model.selected,names(i$analysis$model.selected$data.param))$center2LLLow
             )
           }),
           center=sapply(loci,function(i){
             ifelse(
              is.null(hzar.get.ML.cline(i$analysis$model.selected)$param.all$center),
              NA,
              hzar.get.ML.cline(i$analysis$model.selected)$param.all$center
             )
           }),
           center2LLHigh=sapply(loci,function(i){
             ifelse(
               is.null(hzar.get.ML.cline(i$analysis$model.selected)$param.all$center),
               NA,
               hzar.getLLCutParam(i$analysis$model.selected,names(i$analysis$model.selected$data.param))$center2LLHigh
             )
           }),
           width2LLLow=sapply(loci,function(i){
             ifelse(
               is.null(hzar.get.ML.cline(i$analysis$model.selected)$param.all$width),
               NA,
               hzar.getLLCutParam(i$analysis$model.selected,names(i$analysis$model.selected$data.param))$width2LLLow
             )
           }),
           width=sapply(loci,function(i){
             ifelse(
               is.null(hzar.get.ML.cline(i$analysis$model.selected)$param.all$width),
               NA,
               hzar.get.ML.cline(i$analysis$model.selected)$param.all$width
             )
           }),
           width2LLHigh=sapply(loci,function(i){
             ifelse(
               is.null(hzar.get.ML.cline(i$analysis$model.selected)$param.all$width),
               NA,
               hzar.getLLCutParam(i$analysis$model.selected,names(i$analysis$model.selected$data.param))$width2LLHigh
             )
           }),
           pMin2LLLow=sapply(loci,function(i){
             ifelse(
               is.null(hzar.get.ML.cline(i$analysis$model.selected)$param.all$pMin),
               NA,
               hzar.getLLCutParam(i$analysis$model.selected,names(i$analysis$model.selected$data.param))$pMin2LLLow
             )
           }),
           pMin=sapply(loci,function(i){
             ifelse(
               is.null(hzar.get.ML.cline(i$analysis$model.selected)$param.all$pMin),
               NA,
               hzar.get.ML.cline(i$analysis$model.selected)$param.all$pMin
             )
           }),
           pMin2LLHigh=sapply(loci,function(i){
             ifelse(
               is.null(hzar.get.ML.cline(i$analysis$model.selected)$param.all$pMin),
               NA,
               hzar.getLLCutParam(i$analysis$model.selected,names(i$analysis$model.selected$data.param))$pMin2LLHigh
             )
           }),
           pMax2LLLow=sapply(loci,function(i){
             ifelse(
               is.null(hzar.get.ML.cline(i$analysis$model.selected)$param.all$pMax),
               NA,
               hzar.getLLCutParam(i$analysis$model.selected,names(i$analysis$model.selected$data.param))$pMax2LLLow
             )
           }),
           pMax=sapply(loci,function(i){
             ifelse(
               is.null(hzar.get.ML.cline(i$analysis$model.selected)$param.all$pMax),
               NA,
               hzar.get.ML.cline(i$analysis$model.selected)$param.all$pMax
             )
           }),
           pMax2LLHigh=sapply(loci,function(i){
             ifelse(
               is.null(hzar.get.ML.cline(i$analysis$model.selected)$param.all$pMax),
               NA,
               hzar.getLLCutParam(i$analysis$model.selected,names(i$analysis$model.selected$data.param))$pMax2LLHigh
             )
           }),
           pVal=sapply(loci,function(i){
             ifelse(
               is.null(hzar.get.ML.cline(i$analysis$model.selected)$param.all$pVal),
               NA,
               hzar.get.ML.cline(i$analysis$model.selected)$param.all$pVal
             )
           }),
           row.names=names(loci),
           stringsAsFactors=FALSE
           )
saveRDS(loci,"loci.RDS")
write.csv(cline.data,"cline.data.csv")

#######
#Make plots of all SNP clines
#######


## RED CENTER

### New, Loop through loci to make the overlapping cline plots
cline.data$center2LLHigh-cline.data$center2LLLow -> cline.data$centerCI
top25p_center <- cline.data %>%
  filter(!is.na(centerCI)) %>%
  arrange(centerCI) %>%
  head(floor(94*.25)) %>%
  use_series(X)

pdf("all.loci.confident.center.red.pdf",6,6)
hzar.plot.cline(loci[[1]]$analysis$model.selected,main="",type="n",ylab="SNP frequency",xlab="Pacific coastline (km)")
for (i in names(loci)){
  if (loci[[i]]$analysis$model.name == "nullModel") {next}
  cline_red_RGB <- 0
  if (i %in% top25p_center) {cline_red_RGB <- 1}
  # Could put IF statement here that only plots narrowest CI clines
  # Check if pMin < pMax
  #current_model <- hzar.get.ML.cline(loci[[i]]$analysis$model.selected)
  #if (current_model$param.all$pMax > current_model$param.all$pMin){
  #  current_model$clineFunc <- function(x){
  #hzar.plot.cline(current_model,col=alpha(rgb(cline_red_RBG,0,0), 0.5),pch=NA,add=T) #Just set pch=NA to not plot points, and add=T to plot multiple on one graph. Lines do successfully get darker when they overlap. Some clines go up and others go down.
  hzar.plot.cline(loci[[i]]$analysis$model.selected,col=alpha(rgb(cline_red_RGB,0,0), 0.5),pch=NA,add=T) #Just set pch=NA to not plot points, and add=T to plot multiple on one graph. Lines do successfully get darker when they overlap. Some clines go up and others go down.
}
dev.off()



########
## RED WIDTH


### New, Loop through loci to make the overlapping cline plots
cline.data$width2LLHigh-cline.data$width2LLLow -> cline.data$widthCI
top25p_width <- cline.data %>%
  filter(!is.na(widthCI)) %>%
  arrange(widthCI) %>%
  head(floor(94*.25)) %>%
  use_series(X)

pdf("all.loci.confident.width.red.pdf",6,6)
hzar.plot.cline(loci[[1]]$analysis$model.selected,main="",type="n",ylab="SNP frequency",xlab="Pacific coastline (km)")
for (i in names(loci)){
  if (loci[[i]]$analysis$model.name == "nullModel") {next}
  cline_red_RGB <- 0
  if (i %in% top25p_width) {cline_red_RGB <- 1}
  # Could put IF statement here that only plots narrowest CI clines
  # Check if pMin < pMax
  #current_model <- hzar.get.ML.cline(loci[[i]]$analysis$model.selected)
  #if (current_model$param.all$pMax > current_model$param.all$pMin){
  #  current_model$clineFunc <- function(x){
  #hzar.plot.cline(current_model,col=alpha(rgb(cline_red_RBG,0,0), 0.5),pch=NA,add=T) #Just set pch=NA to not plot points, and add=T to plot multiple on one graph. Lines do successfully get darker when they overlap. Some clines go up and others go down.
  hzar.plot.cline(loci[[i]]$analysis$model.selected,col=alpha(rgb(cline_red_RGB,0,0), 0.5),pch=NA,add=T) #Just set pch=NA to not plot points, and add=T to plot multiple on one graph. Lines do successfully get darker when they overlap. Some clines go up and others go down.
}
dev.off()


#########

##ALL BLACK

pdf("all.loci.black.pdf",6,6)
hzar.plot.cline(loci[[1]]$analysis$model.selected,main="",type="n",ylab="SNP frequency",xlab="Pacific coastline distance (km)",xlim=c(0,4000))
for (i in names(loci)){
  if (loci[[i]]$analysis$model.name == "nullModel") {next}
  hzar.plot.cline(loci[[i]]$analysis$model.selected,col=alpha(rgb(0,0,0), 0.5),pch=NA,add=T) #Just set pch=NA to not plot points, and add=T to plot multiple on one graph. Lines do successfully get darker when they overlap. Some clines go up and others go down.
}
abline(v=c(1843,2674,2916,3386),lty=3,col='darkgray')
mtext(at=c(850,2259,2795,3151,3750),text=c('AK','BC','A','OR','CA'),side=1,line=-1,cex=0.8)
mtext(at=c(2795),text=c('W'),side=1,line=-1.7,cex=0.8)
mtext(at=0,text="North",side=1,line=2)
mtext(at=4000,text="South",side=1,line=2)
dev.off()
########
#additional stuff



## Plot the maximum likelihood cline for the selected model on RStudio graphical device
#for (i in names(loci)){
#  hzar.plot.cline(loci[[i]]$analysis$model.selected,main=loci[[i]]$ref_allele,col=alpha(rgb(0,0,0), 0.5),pch=NA,type="n") #Just set pch=NA to not plot points, and add=T to plot multiple on one graph. Lines do successfully get darker when they overlap. Some clines go up and others go down.
#}



## Plot the 95% credible cline region for the selected model, skip if nullModel
#  Takes pretty long...
#for (i in names(loci)){
#  if (loci[[i]]$analysis$model.name != "nullModel"){
#pdf(paste0("credplots/",loci[[i]]$ref_allele,".pdf"),6,6)
#hzar.plot.fzCline(loci[[i]]$analysis$model.selected,main=loci[[i]]$ref_allele)
#dev.off()
#    }
#}


### Merge cline SNPs with chromosomal mapping data

cline.data <- read.csv("cline.data.csv",stringsAsFactors=F)
SNP_chr <- readRDS("SNP_chr.RDS")
"SNP" -> names(cline.data)[which(names(cline.data)=="X")]
cline.data <- merge(cline.data,SNP_chr,by="SNP",all.x=T)
cline.yes <- cline.data %>%
  filter(model != "nullModel") %>%
  use_series(CHR) %>%
  table %>%
  data.frame %>%
  setNames(nm=c("CHR","Freq"))

cline.all <- cline.data %>%
  use_series(CHR) %>%
  table %>%
  data.frame %>%
  setNames(nm=c("CHR","Freq"))

df2 <- merge(cline.all,cline.yes,by="CHR",all.x=T) %>%
  setNames(nm=c("CHR","Count.all","Count.yes"))

0 -> df2$Count.yes[which(is.na(df2$Count.yes))]

df2$Freq.all <- df2$Count.all/sum(df2$Count.all)
df2$Freq.yes <- df2$Count.yes/sum(df2$Count.yes)

pdf("cline_vs_all.pdf",10,5)
barplot(t(df2[,c("Freq.all","Freq.yes")]),beside=T,names.arg=df2$CHR,las=2,ylab="Frequency") ## (gray=cline,black=all)
dev.off()

chisq.test( df2[,c("Count.all","Count.yes")],
            simulate.p.value=T) ## Not significant

cline.data %>%
  arrange(CHR,BP) %>%
  write.csv("sorted.csv")

## Get some stats for the SNP-based clines

cline.stats <- cline.data %>%
  filter(model != "nullModel")

nrow(cline.stats)

cline.stats$center2LLLow %>% median
cline.stats$center %>% median
cline.stats$center2LLHigh %>% median

cline.stats$width2LLLow %>% median
cline.stats$width %>% median
cline.stats$width2LLHigh %>% median

cline.stats %>% nrow