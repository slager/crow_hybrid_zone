library(ape)
library(magrittr)

nex <- read.nexus.data("allcrowsND2_260.nex")
length(nex) # 260
nex$Corvus_corone <- NULL      ## Remove Outgroup
nex$gws4003_Pierce_Co <- NULL  ## Remove Partial sequence
length(nex) # 258

haps <- read.csv("popart_haplotype259.csv",header=T,stringsAsFactors=F)
haps <- haps[which(haps$X!="Corvus_corone"),]
haps[,"A"] <- as.logical(haps[,"A"])
haps[,"N"] <- as.logical(haps[,"N"])

## Check that they include same names

all(names(nex) %in% haps$X) # TRUE
all(haps$X %in% names(nex)) # TRUE

## Fix the order
haps[match(names(nex), haps$X),] -> haps

## Check the order
all(haps$X == names(nex))   ## TRUE; N=258

## Raw, pairwise, between American and Northwestern
dm <- dist.dna(as.DNAbin(nex),model="raw",as.matrix=T)
dm <- dm[haps[haps$A,"X"],haps[haps$N,"X"]]

nrow(dm) #A
ncol(dm) #N

mean(dm)
range(dm)

median(dm*1041)
range(dm*1041)

### mtDNA Divergece Dating using Lerner et al.
#.029 = subst/site/1e6 yr
#.029 = mean(dm)/(yr/1e6)
# 381 kya (461-335 kya)

mean(dm)/.029*1e6 #estimate
mean(dm)/.024*1e6 #lower HPD
mean(dm)/.033*1e6 #high  HPD

## Identify the fixed differences

haps[which(haps$A==TRUE),'X'] -> amcr.names
haps[which(haps$N==TRUE),'X'] -> nocr.names
nex[amcr.names] -> amcr.nex
nex[nocr.names] -> nocr.nex

sapply(1:length(nex[[1]]),function(x){
  sapply(amcr.nex,function(i){i[[x]]}) -> amcr.nuc
  sapply(nocr.nex,function(i){i[[x]]}) -> nocr.nuc
  if (  (length(unique(amcr.nuc))==1) & (length(unique(nocr.nuc))==1) & (!identical(unique(amcr.nuc),unique(nocr.nuc))) ) return(TRUE) else return(FALSE)
}) %>% unlist %>% which

# ## Make de-duped nexus file
# nex[c(nocr.names,amcr.names)] -> ordered_nex
# ordered_nex[which(!duplicated(ordered_nex))] -> ordered_deduped_nex
# write.nexus.data(ordered_deduped_nex,file="ordered_deduped.nex",interleaved=F)

## Chi-squared test for Vancouver Island vs. Adjacent Mainland
read.csv('pop_haps.csv',stringsAsFactors=F,header=T) -> ph
ph[ph$pop %in% c("nvi","vic","cbc","yvr"),]
matrix(c(25,9,9,26),ncol=2,byrow=T) %>% chisq.test(simulate.p.value=T)
