ref=/media/nfs2/wbs/wbs1/working_dir/genomes/Quercus_mongolica2_update/QuMo.fa
mkdir -p dotplot

minimap2 -t 48 -cx asm10 ${ref} result/ragtag.scaffold.fasta > dotplot/QuMo_$1.paf
