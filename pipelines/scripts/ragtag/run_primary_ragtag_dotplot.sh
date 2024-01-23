ref1=/media/nfs2/wbs/wbs3/working_dir/genomes/SDTS/SDTS.h1.chr.fasta
ref2=/media/nfs2/wbs/wbs3/working_dir/genomes/SDTS/SDTS.h2.chr.fasta
qlen=5e+06

/media/nfs2/wbs/wbs1/bin/paf2dotplot/paf2dotplot.r -s -q ${qlen} \
-i Chr1a,Chr2a,Chr3a,Chr4a,Chr5a,Chr6a,Chr7a,Chr8a,Chr9a,Chr10a,Chr11a,Chr12a \
-o dotplot/SDTS.h1_$1.p_ref1 dotplot/SDTS.h1_$1.p_ref1.paf

/media/nfs2/wbs/wbs1/bin/paf2dotplot/paf2dotplot.r -s -q ${qlen} \
-i Chr1b,Chr2b,Chr3b,Chr4b,Chr5b,Chr6b,Chr7b,Chr8b,Chr9b,Chr10b,Chr11b,Chr12b \
-o dotplot/SDTS.h2_$1.p_ref2 dotplot/SDTS.h2_$1.p_ref2.paf
