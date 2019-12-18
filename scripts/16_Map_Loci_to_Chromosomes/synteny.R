library(magrittr)
library(plyr)
library(dplyr)
library(data.table)
library(tidyr)

## Data frame of all stacks loci
stacks <- fread("cb62_batch_1.sumstats.tsv") %>%
  as.tbl %>%
  select(LOCUS=`Locus ID`,SCAFFOLD=Chr,BP,COL=Col)
stacks$SNP <- paste0(stacks$LOCUS,"_",stacks$COL)
nrow(stacks) # 9591

## Get and concatenate the Coords files produced by Mummer (SKIP AHEAD to readRDS, or uncomment the following to re-run)
# Get only the coords files
# coords.files <- list.files()[which(substr(list.files(),nchar(list.files())-9,nchar(list.files())) == ".fa.coords")]
# Make empty list same length as chromosomes
# coords.list <- vector("list", length(coords.files))
# Make list of coords data frames
# for (i in seq_along(coords.list)){
# coords.list[[i]] <- fread(coords.files[i],skip=5) %>%
  # as.tbl %>%
  # separate(V1,c("S1","E1")) %>%
  # separate(V2,c("S2","E2")) %>%
  # separate(V3,c("LEN1","LEN2")) %>%
  # rename(IDY=V4) %>%
  # separate(V5,c("CHR","SCAFFOLD"),sep="\t")
# }
# Combine list of data frames into one df
# bind_rows(coords.list) -> coords
#Save/Read to RDS
#saveRDS(coords,"coords.RDS")
coords <- readRDS("coords.RDS")   ### READ RDS SHORTCUT ###
coords <- coords %>% rename(SCAFFOLD = CONTIG)
coords$S1 <- as.integer(coords$S1)
coords$E1 <- as.integer(coords$E1)
coords$S2 <- as.integer(coords$S2)
coords$E2 <- as.integer(coords$E2)
coords$LEN1 <- as.integer(coords$LEN1)
coords$LEN2 <- as.integer(coords$LEN2)

## Get locus names from 62 crow Stacks alignment
SNP <- readLines("Cb.s62.r0.05.maf0.02.het0.5.str.tsv",n=2) %>%
  strsplit(split="\t") %>%
  unlist %>%
  `[`(c(-1,-2)) ## Customized for this particular Structure file, above

all(SNP %in% stacks$SNP) ## Yay
# 
# loci <- readLines("Cb.s62.r0.05.maf0.02.het0.5.str.tsv.unlinked",n=1) %>%
#   strsplit(split="\t") %>%
#   unlist %>%
#   `[`(-1) %>%
#   sub("_.+","",.) %>%
#   as.integer

#extra code below for using locus names
##%>%
##  sub("_.+","",.) %>%
##  as.integer
## Also consider this for getting list of all SNPs not just unlinked
#readLines("Cb.s62.r0.05.maf0.02.het0.5.str.tsv",n=2)

SNP_df <- data.frame(SNP=SNP)
SNP_scaffold <- merge(SNP_df,stacks,by='SNP',all.x=T)
nrow(SNP_scaffold) # 9563 SNPs from Structure file
length(unique(SNP_scaffold$LOCUS)) # 7292 loci from Structure file now have scaffolds
length(unique(SNP_scaffold$SCAFFOLD)) #SNPs mapped to 613 different AMCR scaffolds

read.csv("Fst_loci_parentals_ForR.csv",stringsAsFactors=F) -> fst
highfst <- fst[which(fst$Fst > 0.6),"locus"] ## "locus" in this file is actually the full SNP ID
# the following used to be after the above line
## %>% sub("_.+","",.) %>% as.numeric

data.frame(SNP=highfst) -> df_highfst
merge(df_highfst,stacks,by='SNP') -> df_highfst

## Get number of unique chormosomes mapped per scaffold
coords1000_summary <- coords %>%
  filter(LEN1 > 1000 & LEN2 > 1000) %>%
  select(CHR,SCAFFOLD) %>%
  group_by(SCAFFOLD) %>%
  summarize(n(),n_distinct(CHR)) %>%
  as.data.frame
# 
# ## Check example individual scaffold for results of filtering (SKIP)
# coords %>%
#   filter(LEN1 > 1000 & LEN2 > 1000) %>%
#   filter(SCAFFOLD=="NW_008238404.1") %>%
#   arrange(IDY) %>%
#   as.data.frame


## Only use coords when LEN1 > 1000 and LEN2 > 1000
coords1000 <- coords %>% filter(LEN1 > 1000 & LEN2 > 1000) 

#Create blank data frame
scaffold_map <- data.frame(SCAFFOLD=unique(coords1000$SCAFFOLD),CHR=NA,PCT=NA)
#Fill data frame
for (scaffold in unique(coords1000$SCAFFOLD))
     {
  x <- filter(coords1000,SCAFFOLD==scaffold) %>%
    group_by(CHR) %>%
    summarize(sumLEN2=sum(LEN2)) %>%
    as.data.frame
  nchr <- nrow(x)
  chr <- x[which(x$sumLEN2 == max(x$sumLEN2)),'CHR']
  pct <- x[which(x$CHR==chr),'sumLEN2']/sum(x$sumLEN2)
  scaffold_map[scaffold_map$SCAFFOLD==scaffold,'CHR'] <- chr
  scaffold_map[scaffold_map$SCAFFOLD==scaffold,'PCT'] <- pct
}
# Omit matches < 75% 
NA -> scaffold_map[which(scaffold_map$PCT < 0.75),'CHR']

### Get SNP mapping summary stats for crows
SNP_chr <- merge(SNP_scaffold,scaffold_map,by='SCAFFOLD',all.x=T)
SNP_chr$LOCUS %>% unique %>% length # 7,292 loci in Structure file
SNP_chr %>% filter(!is.na(CHR)) %>% use_series(LOCUS) %>% unique %>% length #6,494 loci mapped to CHR via scaffold
SNP_chr %>% filter(!is.na(CHR)) %>% use_series(SCAFFOLD) %>% unique %>% length #327 scaffolds

### Get proportion of mapped SNPs on each chromosome
SNP_chr <- SNP_chr %>% group_by(LOCUS) %>% filter(row_number()==1) %>% ungroup ## Only use 1 SNP per locus
table(SNP_chr$CHR)/sum(table(SNP_chr$CHR)) -> p.crow
p.crow <- as.data.frame(p.crow)
names(p.crow) <- c("CHR","P.CROW")

SNP_chr %>% filter(!is.na(CHR)) %>%
  group_by(LOCUS) %>%
  filter(COL == min(COL)) %>%
  ungroup %>%
  filter(CHR != "chrM") %>%
  nrow # 6492 unlinked SNPs not on MtDNA
SNP_chr %>% filter(!is.na(CHR)) %>%
  filter(CHR == "chrZ") %>%
  group_by(LOCUS) %>%
  filter(COL == min(COL)) %>%
  ungroup %>%
  nrow # 206 unlinked SNPs on Z
206/6492 # 3.17 percent unlinked Z chromosomes SNPs overall
  
#Save as RDS file for SNP-based clines
saveRDS(SNP_chr,"SNP_chr.RDS")

merge(df_highfst,scaffold_map,by='SCAFFOLD') %>%
  filter(CHR == "chrZ") %>%
  nrow
2/35 # Z chromosome SNPs in 35 "most informative"
matrix(c(206,6492,2,35),nrow=2) %>% chisq.test(correct=T,simulate.p.value=T)

highfst_chr <- merge(df_highfst,scaffold_map,by='SCAFFOLD') # To get chromosomes for the ancestry-informative SNPs
# Only 2/35, or 5.7%, of SNPs on the Z. Similar to proportion of Z in Zf genome (7.1%)
write.csv(highfst_chr,"highfst_chr.csv")

## Get Zebra Finch chromosome info
read.csv("ZebraFinch.csv",stringsAsFactors=F) -> zf
zf$SIZE/sum(zf$SIZE) -> zf$P.ZF

#Merge ZF chromosome info with crowRAD chromosome info
SNP_chr <- merge(zf,p.crow,by='CHR',all=T)
SNP_chr$P.ZF_ovr_P.CROW <- SNP_chr$P.ZF/SNP_chr$P.CROW
rownames(SNP_chr) <- SNP_chr$CHR

## Plot 
pdf("barplot.pdf")
foo2 <- SNP_chr[,c("P.ZF","P.CROW")] %>% as.matrix
barplot(height=t(foo2),beside=T,las=2,ylab="Proportion")
dev.off() ## Black = Zebra Finch genome proportion, White = CrowRAD loci proportion

# pdf("synteny.pdf")
# plot(table(foo$CHR),las=2,ylab="# of loci")
# dev.off()

lm(P.ZF_ovr_P.CROW ~ GC,data=SNP_chr) %>% summary
pdf("chromosome_gc.pdf")
plot(P.ZF_ovr_P.CROW ~ GC,data=SNP_chr,ylab="Proportion of Zebra Finch genome / proportion of crow ddRAD loci",xlab="GC content of chromosome")
## "Proportion of Zebra Finch genome / proportion of crow ddRAD loci"
dev.off()


lm(P.ZF_ovr_P.CROW ~ GENE,data=SNP_chr) %>% summary
plot(P.ZF_ovr_P.CROW ~ GENE,data=SNP_chr)

lm(GENE ~ GC, data=SNP_chr) %>% summary
plot(GENE ~ GC, data=SNP_chr)