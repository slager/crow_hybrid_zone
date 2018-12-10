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
  if (use_locus_names==TRUE){sitelabels -> colnames(s)}
  ## s is now a data.frame
  ## with {#samples*2} rows and {#sites} columns.
  ## Data = various integers
  ## Column names of s are the site labels
  cat(ncol(s),"SNPs\n")
  
  #assign("s", s, envir = .GlobalEnv)
  return(s)
} # End of function parseStructureFile


retainParsimonyInformative <- function(s){
  ### Detect parsimony-informative sites in s
  
  # Create L, a list of vectors for indexing:
  # length(L) = #samples
  # L[[1]] = c(1,2)
  # L[[nSamples]] = c(2*nSamples-1,2*nSamples)
  lapply(seq(1,nrow(s)-1,2),function(x){x:(x+1)}) -> L
  
  #For each site (each column of s):
  apply(s,2,function(a){
    #Sort the diploid unlinked SNP data so that (1,2) is same as (2,1)
    lapply(1:length(L),function(x){sort(a[L[[x]]])}) -> g
    #Remove samples with missing data before detecting pi sites
    g[which(sapply(g,function(x){!all(x==c(missingdatainteger,missingdatainteger))}))] -> g
    #Coerce to character & paste (result = e.g. "11", "12", "22", "13", ...)
    lapply(g,function(x){paste(as.character(x),sep="",collapse="")}) -> g
    unlist(g) -> g  # Convert to vector
    #Check if parsimony informative
    #https://en.wikipedia.org/wiki/Informative_site
    rle(sort(g))$lengths -> sitepatterns
    if (length(sitepatterns[which(sitepatterns >= 2)]) >= 2) TRUE else FALSE #This is value returned by apply() for each site
  }) -> p # A vector of length nsites, indicates if site parsimony-informative
  
  ## Get parsimony-informative data.frame
  s[,which(p)] -> pi
  
  cat(ncol(pi),"parsimony-informative SNPs retained\n")
  #assign("s", pi, envir = .GlobalEnv)
  return(pi)
} # End function retainParsimonyInformative



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

retainRandomUnlinkedSNP <- function(s){
  if (use_locus_names==FALSE){cat("Error: RandomSNP currently requires STACKS-style locus names\n");stop}
  ## Get 1 random unlinked SNP per locus
  sapply(strsplit(colnames(s),"_"),function(x){x[1]}) -> loci
  unique(loci) -> uloci
  cat(ncol(s),"SNPs at",length(uloci),"loci\n")
  rep(as.character(NA),length(uloci)) -> u
  for (i in 1:length(uloci)){
    colnames(s)[which(loci==uloci[i])] -> current_locus_SNPs
    sample(current_locus_SNPs,1)-> u[i]
  }
  s[,u] -> unlinked
  #assign("s", unlinked, envir = .GlobalEnv)
  cat(ncol(unlinked),"random unlinked SNPs retained\n")
  return(unlinked)
} # End of retainRandomUnlinkedSNP

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

# 
# ### OBSOLETE CODE BELOW
# 
# ###  Where the action happens ###
# 
# ## Optional single file action
# #filename <- "str_files/Cb.s64.r0.05.maf0.02.het0.5.str.tsv"
# 
# pops2keep <- c("ca","cbc","ghc","hmr","jun","kit","nbc","neah","nvi","sea","vic","yvr")
# 
# ## Optional "for" loop action
# filenames <- c("Cb.s62.r0.05.maf0.02.het0.5.str.tsv","Cb.s64.r0.05.maf0.02.het0.5.str.tsv")
# directory <- "str_files/" ## Optional
# filenames <- paste0(directory,filenames) ## Optional
# for (filename in filenames){
#   parseStructureFile() #Reads 'filename', writes 's'
#   retainParsimonyInformative()
#   retainRandomUnlinkedSNP() #Reads 's', writes 's'
#   writeOutStructure() #Reads 's', writes 'filename'+'suffix'
# }
# 
# #### OBSOLETE CODE ####
# 
# filename <- "str_files/Cb.s62.r0.05.maf0.02.het0.5.str.tsv"
# parseStructureFile()
# pops2keep <- c("ca","cbc","ghc","hmr","jun","kit","nbc","neah","nvi","sea","vic","yvr")
# selectSamples(pops2keep=pops2keep)
# retainFullCoverage()
# retainRandomUnlinkedSNP()
# writeOutStructure(suffix=".coastal.unlinked.nomd")
# 
# 
# 
# 
# 
# ## Summary outputs to console
# #ncol(s)          # Total number of sites
# #length(unique(sapply(strsplit(colnames(s),"_"),function(x){x[1]}))) # Total number of loci
# #ncol(pi)         # Total number of parsimony informative sites
# #ncol(upi)        # Number of unlinked parsimony informative sites
# #ncol(pi)/ncol(s) # Proportion of sites parsimony informative
# 
# ## Summary outputs to console OLD
# #ncol(s)          # Number of sites
# #ncol(pi)         # Number of parsimony informative sites
# #ncol(pi)/ncol(s) # Proportion of sites parsimony informative
# 
# ## Write New Structure file PI
# # Write site names header (Optional)
# #write(paste(c("",names(pi)),collapse="\t",sep=""),file=paste0(filename,".pi"))
# # Write Structure data
# #write.table(pi,file=paste0(filename,".pi"),col.names=FALSE,row.names=samplelabels,append=TRUE,quote=FALSE,sep="\t")
# 
# 
# #}  #Bracket for entire "for" loop
# 
# 
# ## Older remnant code (FYI) for initial parse of pyRAD STRUCTURE files
# 
# #fread("c88m48p9.str",data.table=F,verbose=T) -> s
# #paste0(s$V1,c("a","b")) -> rownames(s)
# #s[,-c(1:6)] -> s
# #str(s)