library(ape)

#read.nexus.data("allcrowsND2.nex") -> nexus
read.nexus.data("allcrowsND2_260.nex") -> nexus
names(nexus) -> nexus.labels
#write.csv(nexus.labels,"nexus.labels.csv")

as.DNAbin(nexus) -> alignment
dist.dna(alignment) -> distances
bionj(distances) -> nj

#write.nexus(nj,file="nexus.translated",translate=TRUE)
#write.nexus(nj,file="nexus.untranslated",translate=FALSE)
#write.tree(nj,"newick")

root(nj,"Corvus_corone") -> rootednj
makeNodeLabel(rootednj, method = "number", prefix = "N") -> rootednj

# Plot with node labels
pdf("nj_tree.pdf",20,100)
plot.phylo(rootednj,show.node.label=TRUE)
dev.off()

extract.clade(rootednj,"N127")$tip.label -> nw # 95
extract.clade(rootednj,"N123")$tip.label -> am # 164
data.frame(tip.label=c(nw,am),haplotype=c(rep("N",length(nw)),rep("A",length(am))),stringsAsFactors=FALSE) -> info
write.csv(info,"haplotype.csv")
