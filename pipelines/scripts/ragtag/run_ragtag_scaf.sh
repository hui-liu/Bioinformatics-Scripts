source ~/.bash_conda
conda activate ragtag

ref1=/media/nfs2/wbs/wbs3/working_dir/genomes/SDTS/SDTS.h1.chr.fasta
ref2=/media/nfs2/wbs/wbs3/working_dir/genomes/SDTS/SDTS.h2.chr.fasta

# hap1 scaffolding
ragtag.py scaffold -t 48 -o hap1_1 ${ref1} $1.h1.fasta
ragtag.py scaffold -t 48 -o hap1_2 ${ref2} $1.h1.fasta

# hap2 scaffolding
ragtag.py scaffold -t 48 -o hap2_1 ${ref1} $1.h2.fasta
ragtag.py scaffold -t 48 -o hap2_2 ${ref2} $1.h2.fasta
