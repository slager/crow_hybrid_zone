#!/bin/bash

# Filename list of demultiplexed reads without extensions
files=$(ls /home/burke/Documents/dls/crowRAD/stacks/fq.gz/ | sed -e 's/\..*//')

## Align to Corvus cornix genome

for file in $files 
do
	gsnap -D /home/burke/Documents/dls/crowRAD/stacks/gmapdb/ \
		-d Corvus_cornix \
		--gunzip -t 32 -A sam \
		--min-coverage=0.90 --max-mismatches=3 \
		--indel-penalty=2   \
		/home/burke/Documents/dls/crowRAD/stacks/fq.gz/${file}.fq.gz \
			> /home/burke/Documents/dls/crowRAD/stacks/ref_align/Corvus_cornix2/sam/${file}.sam 
	samtools view -b -S \
		-o /home/burke/Documents/dls/crowRAD/stacks/ref_align/Corvus_cornix2/bam/${file}.bam \
		/home/burke/Documents/dls/crowRAD/stacks/ref_align/Corvus_cornix2/sam/${file}.sam
done

## Align to American Crow genome

for file in $files 
do
	gsnap -D /home/burke/Documents/dls/crowRAD/stacks/gmapdb/ \
		-d Corvus_brachyrhynchos \
		--gunzip -t 32 -A sam \
		--min-coverage=0.90 --max-mismatches=3 \
		--indel-penalty=2   \
		/home/burke/Documents/dls/crowRAD/stacks/fq.gz/${file}.fq.gz \
			> /home/burke/Documents/dls/crowRAD/stacks/ref_align/Corvus_brachyrhynchos2/sam/${file}.sam 
	samtools view -b -S \
		-o /home/burke/Documents/dls/crowRAD/stacks/ref_align/Corvus_brachyrhynchos2/bam/${file}.bam \
		/home/burke/Documents/dls/crowRAD/stacks/ref_align/Corvus_brachyrhynchos2/sam/${file}.sam
done
