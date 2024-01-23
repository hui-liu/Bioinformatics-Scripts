ref=/media/nfs2/wbs/wbs1/working_dir/genomes/Quercus_mongolica2_update/QuMo.fa
qlen=5e+06

/media/nfs2/wbs/wbs1/bin/paf2dotplot/paf2dotplot.r -s -q ${qlen} \
-i Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12 \
-o dotplot/QuMo_$1 dotplot/QuMo_$1.paf
