#!/bin/bash
# Scaffold alignment pipeline for avian whole genomes, modified script from CJ Battey
# run time is ~12hrs. 

cd ~/Documents/dls/synteny/synteny
sudo apt install genometools
conda install mummer

#download repeat-masked zebra finch genome and split fasta entries. 
curl "http://hgdownload.cse.ucsc.edu/goldenPath/taeGut2/bigZips/taeGut2.fa.masked.gz" -o "tgut2.fa.gz"
gunzip tgut2.fa.gz
mkdir tgut2_split
gt splitfasta -splitdesc tgut2_split/ tgut2.fa
#rm tgut2.*  
cd ./tgut2_split/ #remove short contigs, leave the chromosome-level scaffolds. (?) 
find . -name "*_*" -delete
cd ~/Documents/dls/synteny/synteny

#download americna crow genome.
curl "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/691/975/GCF_000691975.1_ASM69197v1/GCF_000691975.1_ASM69197v1_genomic.fna.gz" -o "Cb.fa.gz"
gunzip Cb.fa.gz

#output directories
mkdir mummer_out
mkdir ./mummer_out/coords
mkdir ./mummer_out/alignments

#Run mummer alignment against each chromosome. 
refSeqs=tgut2_split/*

for chr in $refSeqs

	do
	echo "aligning to $chr"
	
	nucmer $chr ./Cb.fa
	
	show-coords ./out.delta > ./mummer_out/coords/${chr##*/}.coords #Check regex if changing paths.
	
	mv out.delta ./mummer_out/alignments/${chr##*/}.delta
	
done

