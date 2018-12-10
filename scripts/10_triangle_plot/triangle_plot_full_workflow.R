library(plyr)
library(dplyr)

##
## 1.  In Excel or R, already defined parentals individuals using 98%  CLUMPP cutoffs K=2.
## (Annotated in "parental define from Clumpp 0.98.csv")
##

##
## 2.  parentals_define_structure_functions.R   (just defines functions for subsetting)
##
##

## Define all the functions

parseStructureFile <- function(){
  cat("Parsing",filename,"\n")
  cat("Missing data coded as",missingdatainteger,"\n")
  
  ## Get locus labels vector, if applicable
  if (use_locus_names==TRUE){
    readLines(filename,n=locusnameline)[locusnameline] -> a
    unlist(strsplit(a,split=locusnamesep)) -> sitelabels
    sitelabels[which(sitelabels!="")] -> sitelabels
  }
  
  fread(filename,data.table=F,verbose=F) -> s
  
  ## Get sample labels
  s$V1 -> samplelabels  # Save the original sample labels
  assign("samplelabels", samplelabels, envir = .GlobalEnv)
  # Relabel samples because duplicate rownames not allowed
  paste0(s$V1,c("a","b")) -> rownames(s) #Rownames useful for troubleshooting
  cat("Assuming diploid data sorted by sample ID\n")
  cat(nrow(s)/2,"samples\n")
  
  ## Remove leading columns so that #columns = #sites
  ## This part is file-specific.
  s[,first_locus_column:ncol(s)] -> s
  
  ## Assign locus labels as the column names
  #if (use_locus_names==TRUE){sitelabels -> colnames(s)}
  ## s is now a data.frame
  ## with {#samples*2} rows and {#sites} columns.
  ## Data = various integers
  ## Column names of s are the site labels
  cat(ncol(s),"SNPs\n")
  
  #assign("s", s, envir = .GlobalEnv)
  return(s)
} # End of function parseStructureFile

retainFullCoverage <- function(s){
  ### Remove SNPs with missing data for any individual
  
  #For each site (each column of s):
  apply(s,2,function(a){
    ## TRUE if does not contain missing data
    !(missingdatainteger %in% a) -> no_missing_data
    (length(unique(a)) > 1) -> not_invariant
    as.logical(no_missing_data*not_invariant)
  }) -> toRetain # A vector of length nsites, indicates if site has no missing data
  
  ## Get data.frame of no missing data SNPs
  s[,which(toRetain)] -> nomd
  
  cat(ncol(nomd),"Full coverage SNPs retained\n")
  #assign("s", nomd, envir = .GlobalEnv)
  return(nomd)
} # End function retainFullCoverage

selectSamples <- function(s){
  cat("Assuming sample names format e.g. 'ca01' 'ca02' from population 'ca'\n")
  gsub("[0-9]","",samplelabels) -> poplabels
  kept <- s[which(poplabels %in% pops2keep),]
  samplelabels <- samplelabels[which(poplabels %in% pops2keep)]
  cat(ncol(kept),"SNPs retained\n")
  cat(nrow(kept)/2,"samples retained from",length(pops2keep),"populations\n")
  assign("samplelabels", samplelabels, envir = .GlobalEnv)
  #assign("s", kept, envir = .GlobalEnv)
  return(kept)
}

writeOutStructure <- function(s){
  ## Write New Structure file
  if (nchar(suffix)==0) {warning("Write terminated: Check writeOutSuffix");stop}
  # Write site names header (Optional)
  if (use_locus_names==TRUE){
    write(paste(c("",names(s)),collapse="\t",sep=""),file=paste0(filename,suffix))
  }
  # Write Structure data
  write.table(s,file=paste0(filename,suffix),col.names=FALSE,row.names=samplelabels,append=use_locus_names,quote=FALSE,sep="\t")
  cat("Wrote to",paste0(filename,suffix),"\n")
} #End of writeOutStructure

##
## 3.  parentals_Subset_Structure.R  (just subsets to the parental individuals, defined by parentalA and parentalN vectors)
## Input is Cb.s62.r0.05.maf0.02.het0.5.tsv.unlinked.str
## Output is Cb.s62.r0.05.maf0.02.het0.5.tsv.unlinked.str.parentals.str)

library(data.table)
library(magrittr)

missingdatainteger <- 0     # 0 stacks, -9 pyRAD
locusnameline <- 2    # Line number of locus names in the infile
locusnamesep <- "\t"  # Delimiter for locus names line
use_locus_names <- TRUE # Read/use/write STACKS-style locus names?
first_locus_column <- 2 #Column no. of 1st SNP data in Structure file
suffix <- ".parentals.str"
pops2keep <- c("ca","cbc","ewa","ghc","hmr","jun","kit","la","mi","nbc","neah","nf","nvi","nynj","sea","vic","yvr")

filenames <- c("Cb.s62.r0.05.maf0.02.het0.5.tsv.unlinked.str")
directory <- "" ## Optional
filename <- paste0(directory,filenames) ## Optional

parseStructureFile() -> s

#Get normal sample labels back (remove a/b)
s$id <- substr(rownames(s),1,nchar(rownames(s))-1)

# Read custom CSV of Clumpp-defined parentals >98% pure
parental_df <- fread("parental define from Clumpp 0.98.csv",data.table=F)

# Sample IDs for parentals only
parentalN <- parental_df %>% filter(N_parental==1) %>% use_series(id)
parentalA <- parental_df %>% filter(A_parental==1) %>% use_series(id)


#Subset structure file to only parentals
s <- filter(s,id %in% c(parentalA,parentalN))

#Save sample IDs for later
s$id -> sid_save
#Delete sample IDs
s <- s %>% select(-id)
#Set 0 to NA
as.matrix(s) -> s
NA -> s[which(s==0)]
as.data.frame(s) -> s
#Number of samples w/ data at a SNP (diploid); max 28 (parentals)
vapply(s,FUN.VALUE=1,function(x){length(na.omit(x))})/2 -> vec
# 28 * .75 = 21
s <- s[,vec/28>=0.75] ## Require 75% of parentals to have data for locus
ncol(s)
#Remove invariant SNPs
vapply(s,FUN.VALUE=TRUE,function(x){
  length(unique(na.omit(x)))>1}) -> is.variable.site
s[,is.variable.site] -> s
ncol(s)
#Set NA to 0
as.matrix(s) -> s
0 -> s[which(is.na(s))]
as.data.frame(s) -> s
# Re-add sample names for rownames
#sid_save -> s$id
#s <- select(s,id,everything())

sid_save -> samplelabels
#writeOutStructure(s)     ## Want this ON if running 1st time



##
## 4.  split_parentals_structure_and_make_Fst.R
## (reads in *parentals.str and CLUMPP 98% purity definitions, write Fst file)
##

library(data.table)
library(magrittr)
library(plyr)
library(dplyr)
library(adegenet)
library(hierfstat)

# Read structure file to data.table
s <- fread("Cb.s62.r0.05.maf0.02.het0.5.tsv.unlinked.str.parentals.str")


# Read big structure file to genind
read.structure("Cb.s62.r0.05.maf0.02.het0.5.tsv.unlinked.str.parentals.str",
               n.ind=28,n.loc=3582,onerowperind=F,
               col.lab=1,col.pop=0,col.others=0,
               row.marknames=1,NA.char="0") -> str


#Subset Genind object
str[c(parentalN,parentalA),] -> str
pop_vector <- c(rep(1,length(parentalN)),rep(2,length(parentalA)))

genind2hierfstat(str,pop=pop_vector) -> str_fstat

# Structure files for parentals only
parentalN_str <- s %>% filter(V1 %in% parentalN)
parentalA_str <- s %>% filter(V1 %in% parentalA)

#basic.stats(str_fstat)$perloc %>% write.csv("Fst_loci_parentals.csv")  ## Actually write this file, if running 1st time
#basic.stats(str_fstat)$perloc$Fst %>% hist(breaks=(-2:20)/20,ylim=c(0,50))



##
## 5.  informative_loci_Subset_Structure.R
## starts with original 62-sample file, because eventually plot all 62 indiv.
## (keeps informative loci, based on file Fst_loci_parentals_ForR.csv)
##

library(data.table)
library(magrittr)

missingdatainteger <- 0     # 0 stacks, -9 pyRAD
locusnameline <- 2    # Line number of locus names in the infile
locusnamesep <- "\t"  # Delimiter for locus names line
use_locus_names <- TRUE # Read/use/write STACKS-style locus names?
first_locus_column <- 2 #Column no. of 1st SNP data in Structure file
suffix <- ".informative.loci.str"
pops2keep <- c("ca","cbc","ewa","ghc","hmr","jun","kit","la","mi","nbc","neah","nf","nvi","nynj","sea","vic","yvr")

filenames <- c("Cb.s62.r0.05.maf0.02.het0.5.tsv.unlinked.str")
directory <- "" ## Optional
filename <- paste0(directory,filenames) ## Optional
#filename <- "/Users/dave/crows/stacks/missing_data/str_files/Cb.s62.r0.05.maf0.02.het0.5.str.tsv"

parseStructureFile() -> s

#Get normal sample labels back (remove a/b)
s$id <- substr(rownames(s),1,nchar(rownames(s))-1)

#Subset structure file to only parentals
#s <- filter(s,id %in% c(parentalA,parentalN))

#Subset structure file to only informative loci
#import informative loci file
fread("Fst_loci_parentals_ForR.csv") -> fst_file
fst_file <- fst_file %>% filter(Fst>0.6) %>% select(locus,Fst)
informative_loci <- fst_file$locus
#informative_loci <- c("13037_27","20857_41","8613_26")

# Prepare admix.gen data frame for input
s <- select(s,c("id",informative_loci)) #informative loci only
admix.gen_str <- s
#admix.gen_str <- filter(s, ! id %in% c(parentalA,parentalN))
#dim(admix.gen_str)

# Subset parental only structure files
str_pN <- s %>% filter(id %in% parentalN) %>% select(-id)
dim(str_pN)
str_pA <- s %>% filter(id %in% parentalA) %>% select(-id)
dim(str_pA)

library(introgress)

## Get data into Introgress
introgress_prepared_data <-
  prepare.data(
    admix.gen=t(admix.gen_str),
    parental1=t(str_pN),
    parental2=t(str_pA),
    loci.data=data.frame(locus.name=colnames(str_pN),type=rep("C",ncol(str_pN))),
    pop.id=FALSE,
    ind.id=TRUE,
    fixed=FALSE,
    sep.columns=TRUE)

est.h(
  introgress.data=introgress_prepared_data,
  loci.data=data.frame(locus.name=colnames(str_pN),type=rep("C",ncol(str_pN))),
  ind.touse=NULL, #use all Indiv if NULL
  fixed=FALSE,
  p1.allele=NULL,
  p2.allele=NULL) -> hi_df

calc.intersp.het(
  introgress.data=introgress_prepared_data) -> het_df

#Package triangle plot
introgress::triangle.plot(
  hi.index=hi_df,
  int.het=het_df,
  pdf=FALSE)

#Plot data frame
data.frame(id=(admix.gen_str$id %>% unique),
           hi_lower=hi_df$lower,
           hi=hi_df$h,
           hi_upper=hi_df$upper,
           het=het_df) %>% arrange(hi) -> df

#df

## Custom colors to match map
'#CC3311' -> red
'#0077BB' -> blue

#My pretty custom triangle plot
pdf("triangle2.pdf",5.5,5.5) #Calling it triangle2 now that it has map theme colors customized
par(mar=c(4.1,4.1,4.1,2.6))
plot(1:10,1:10,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Genomic hybrid index",ylab="Heterozygosity")
segments(0,0,.5,1,lty=5)
segments(.5,1,1,0,lty=5)
palette <- colorRampPalette(c(red,blue))(100)
points(1-df$hi,df$het,pch=21,bg=palette[cut(df$hi,100)]) ## Plotting 1 minus hybrid index to make plot more intuitive
# text(
#   (df %>% filter(id=="nbc04"))$hi,
#   (df %>% filter(id=="nbc04"))$het,
#   labels="nbc04",
#   pos=4,
#   offset=1
#   )
mtext("Northwestern",side=1,line=2,at=0,cex=0.9)
mtext("American",side=1,line=2,at=1,cex=0.9)
dev.off()