#!/bin/bash

# Set Bam directory
B="/home/burke/Documents/dls/crowRAD/stacks/ref_align/Corvus_brachyrhynchos/bam"


# Corvus_brachyrhynchos 64 samples


ref_map.pl \
-o /home/burke/Documents/dls/crowRAD/stacks/output/Corvus_brachyrhynchos.s64.r0.05.maf0.02.het0.5/ \
-T 32 \
-b 1 -S \
-s $B/ca01.bam \
-s $B/ca02.bam \
-s $B/ca03.bam \
-s $B/ca04.bam \
-s $B/cbc01.bam \
-s $B/cbc02.bam \
-s $B/cbc03.bam \
-s $B/cbc04.bam \
-s $B/ewa01.bam \
-s $B/ewa02.bam \
-s $B/ewa03.bam \
-s $B/ewa04.bam \
-s $B/ewa05.bam \
-s $B/ghc01.bam \
-s $B/ghc02.bam \
-s $B/ghc03.bam \
-s $B/ghc04.bam \
-s $B/hmr01.bam \
-s $B/hmr02.bam \
-s $B/hmr03.bam \
-s $B/hmr04.bam \
-s $B/jun01.bam \
-s $B/jun02.bam \
-s $B/jun03.bam \
-s $B/jun04.bam \
-s $B/kit01.bam \
-s $B/kit02.bam \
-s $B/kit03.bam \
-s $B/kit04.bam \
-s $B/la01.bam \
-s $B/la02.bam \
-s $B/mi01.bam \
-s $B/mi02.bam \
-s $B/nbc01.bam \
-s $B/nbc02.bam \
-s $B/nbc03.bam \
-s $B/nbc04.bam \
-s $B/neah01.bam \
-s $B/neah02.bam \
-s $B/neah03.bam \
-s $B/neah04.bam \
-s $B/nf01.bam \
-s $B/nvi01.bam \
-s $B/nvi02.bam \
-s $B/nvi03.bam \
-s $B/nvi04.bam \
-s $B/nynj01.bam \
-s $B/nynj02.bam \
-s $B/nynj03.bam \
-s $B/nynj04.bam \
-s $B/rus01.bam \
-s $B/rus02.bam \
-s $B/sea01.bam \
-s $B/sea02.bam \
-s $B/sea03.bam \
-s $B/sea04.bam \
-s $B/vic01.bam \
-s $B/vic02.bam \
-s $B/vic03.bam \
-s $B/vic04.bam \
-s $B/yvr01.bam \
-s $B/yvr02.bam \
-s $B/yvr03.bam \
-s $B/yvr04.bam \
-X "populations:-e sbfI" -X "populations:-r 0.05" -X "populations:--min_maf 0.02" -X "populations: --max_obs_het 0.5" -X "populations: -k " -X "populations:--ordered_export" -X "populations:--genomic" -X "populations: --structure "  -X "populations: --hzar " -X "populations:--fstats" -X "populations:--fasta" -X "populations:--fasta_strict" -X "populations: --phylip " -X "populations: --phylip_var" -X "populations: --phylip_var_all "


# Corvus_brachyrhynchos 62 samples


ref_map.pl \
-o /home/burke/Documents/dls/crowRAD/stacks/output/Corvus_brachyrhynchos.s62.r0.05.maf0.02.het0.5/ \
-T 32 \
-b 1 -S \
-s $B/ca01.bam \
-s $B/ca02.bam \
-s $B/ca03.bam \
-s $B/ca04.bam \
-s $B/cbc01.bam \
-s $B/cbc02.bam \
-s $B/cbc03.bam \
-s $B/cbc04.bam \
-s $B/ewa01.bam \
-s $B/ewa02.bam \
-s $B/ewa03.bam \
-s $B/ewa04.bam \
-s $B/ewa05.bam \
-s $B/ghc01.bam \
-s $B/ghc02.bam \
-s $B/ghc03.bam \
-s $B/ghc04.bam \
-s $B/hmr01.bam \
-s $B/hmr02.bam \
-s $B/hmr03.bam \
-s $B/hmr04.bam \
-s $B/jun01.bam \
-s $B/jun02.bam \
-s $B/jun03.bam \
-s $B/jun04.bam \
-s $B/kit01.bam \
-s $B/kit02.bam \
-s $B/kit03.bam \
-s $B/kit04.bam \
-s $B/la01.bam \
-s $B/la02.bam \
-s $B/mi01.bam \
-s $B/mi02.bam \
-s $B/nbc01.bam \
-s $B/nbc02.bam \
-s $B/nbc03.bam \
-s $B/nbc04.bam \
-s $B/neah01.bam \
-s $B/neah02.bam \
-s $B/neah03.bam \
-s $B/neah04.bam \
-s $B/nf01.bam \
-s $B/nvi01.bam \
-s $B/nvi02.bam \
-s $B/nvi03.bam \
-s $B/nvi04.bam \
-s $B/nynj01.bam \
-s $B/nynj02.bam \
-s $B/nynj03.bam \
-s $B/nynj04.bam \
-s $B/sea01.bam \
-s $B/sea02.bam \
-s $B/sea03.bam \
-s $B/sea04.bam \
-s $B/vic01.bam \
-s $B/vic02.bam \
-s $B/vic03.bam \
-s $B/vic04.bam \
-s $B/yvr01.bam \
-s $B/yvr02.bam \
-s $B/yvr03.bam \
-s $B/yvr04.bam \
-X "populations:-e sbfI" -X "populations:-r 0.05" -X "populations:--min_maf 0.02" -X "populations: --max_obs_het 0.5" -X "populations: -k " -X "populations:--ordered_export" -X "populations:--genomic" -X "populations: --structure "  -X "populations: --hzar " -X "populations:--fstats" -X "populations:--fasta" -X "populations:--fasta_strict" -X "populations: --phylip " -X "populations: --phylip_var" -X "populations: --phylip_var_all "


