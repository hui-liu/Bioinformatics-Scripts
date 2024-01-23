ref1=/media/nfs2/wbs/wbs3/working_dir/genomes/SDTS/SDTS.h1.chr.fasta
ref2=/media/nfs2/wbs/wbs3/working_dir/genomes/SDTS/SDTS.h2.chr.fasta

minimap2 -t 48 -cx asm10 ${ref1} hap1_1/ragtag.scaffold.fasta > dotplot/SDTS.h1_$1.h1.paf
minimap2 -t 48 -cx asm10 ${ref2} hap1_2/ragtag.scaffold.fasta > dotplot/SDTS.h2_$1.h1.paf
minimap2 -t 48 -cx asm10 ${ref1} hap2_1/ragtag.scaffold.fasta > dotplot/SDTS.h1_$1.h2.paf
minimap2 -t 48 -cx asm10 ${ref2} hap2_2/ragtag.scaffold.fasta > dotplot/SDTS.h2_$1.h2.paf
