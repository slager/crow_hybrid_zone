#setwd("/Users/dave/crows/stacks/missing_data/")
library(data.table)
library(magrittr)

missingdatainteger <- 0     # 0 stacks, -9 pyRAD
locusnameline <- 2    # Line number of locus names in the infile
locusnamesep <- "\t"  # Delimiter for locus names line
use_locus_names <- TRUE # Read/use/write STACKS-style locus names?
first_locus_column <- 3 #Column no. of 1st SNP data in Structure file
suffix <- ".pi.nomd.str2"
pops2keep <- c("ca","cbc","ghc","hmr","jun","kit","nbc","neah","nvi","sea","vic","yvr")

filenames <- c("Cb.s62.r0.05.maf0.02.het0.5.str.tsv","Cb.s64.r0.05.maf0.02.het0.5.str.tsv")
directory <- "str_files/" ## Optional
filenames <- paste0(directory,filenames) ## Optional
#filename <- "/Users/dave/crows/stacks/missing_data/str_files/Cb.s62.r0.05.maf0.02.het0.5.str.tsv"


## Functions available & recommended order
#parseStructureFile()
#selectSamples()
#retainParsimonyInformative()
#retainFullCoverage()
#retainRandomUnlinkedSNP()
#writeOutStructure()

## String desired functions together with the '%>%' pipe

for (filename in filenames){

  parseStructureFile() %>%
  retainParsimonyInformative %>%
  retainRandomUnlinkedSNP %>%
  writeOutStructure()
  
}