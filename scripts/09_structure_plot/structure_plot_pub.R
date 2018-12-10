library(data.table)
library(magrittr)

fread("K2.indfile.outfile.csv",header=F,data.table=F) -> s
fread("haplotypes.csv") -> h

s$V1 -> rownames(s)

s[,-1] -> s ## Remove "indiv" column
c("A","N") -> colnames(s)
str(s)

## Resulting data frame contains N rows, K columns, rownames=sample IDs
#s

rownames(s) %>% sub("ca","nca",.) -> rownames(s)


indexing <- c(18:25,34:37,5:8,43:46,55:62,38:41,26:29,51:54,14:17,1:4,9:13,30:33,47:50,42)

haps <- sapply(rownames(s),
               function(x){h[which(h$rad.id==x),'haplotype.NJ']}
               ) %>% unname %>% unlist
isA <- haps %in% c('A')

s[c(indexing),] -> s
haps <- haps[indexing]
isA <- isA[indexing]

62 -> n
x.txt <- (1:n)-0.5

## Get label midpoint values

pop.numbers <- rownames(s) %>%
  substr(nchar(.),nchar(.))

pop.names <- rownames(s) %>%
  gsub("[012345]","",.) %>%
  unique

pop.values <- pop.names %>%
  sapply(.,function(x){which((rownames(s) %>% gsub("[012345]","",.))==x)})

group.at <-pop.values %>%
  lapply(mean) %>%
  unlist
group.at <- 62.5 - group.at

pop.values <- lapply(pop.values,function(x){62.5-x})
pop.values <- lapply(pop.values,function(x){
  if (length(x)==1) {return(c(x-.25,x+.25))} else {return(c(min(x)-.25,max(x)+.25))}})

## Stacked bar plot grouped by locality
#graph presets
haplotype_x <- 2.15
text_x <- -.05
brace_x <- -.1
brace_width <- .025
#popnum_x <- -.05  # no longer plotting population numbers

'#CC3311' -> red
'#0077BB' -> blue


#png("custom3.png",w=13/6,h=6.5,units="in",res=1000)
pdf(paste0("structure_plot_final.pdf"),13/6,6.5,useDingbats=F)
par(mar=rep(0,4))
barplot(2*as.table(t(s[nrow(s):1,2:1])),xaxt="n",yaxt="n",xaxs="i",yaxs="i",space=0,col=c(blue,red),main='',horiz=T,xlim=c(-0.7,2.3),ylim=c(-0.1,62.1))
#axis(1,at=(1:n)-.5,labels=rownames(s),lwd=0,lwd.ticks=0,las=2) # or lwd.ticks=1
#mtext(at=x.txt,side=4,line=0,text=rev(haps),col=rev(ifelse(isA,'red','blue')))
#mtext(at=x.txt,side=4,line=0,text=".",,padj=-.1,col=rev(ifelse(isA,'red','blue')),cex=3.5)
points(rep(haplotype_x,nrow(s)),x.txt,pch=21,bg=rev(ifelse(isA,red,blue)))
text(x=text_x,y=group.at-.1,labels=pop.names,pos=2,cex=0.8)
for (i in 1:length(pop.values)){
  segments(brace_x,pop.values[[i]][1],brace_x,pop.values[[i]][2],ljoin=1)
  segments(brace_x,unlist(pop.values),brace_x+brace_width,unlist(pop.values),ljoin=1)
}
#text(x=popnum_x,y=x.txt,labels=rev(pop.numbers),cex=0.75)
#mtext(rev(unique(rownames(s))),side=2,
#abline(v=seq(4,44,by=4),lwd=2,col="black") # Lines to separate my localities (line optional)
dev.off()