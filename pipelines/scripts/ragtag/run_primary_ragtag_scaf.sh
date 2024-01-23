source ~/.bash_conda
conda activate ragtag

ref1=/media/nfs2/wbs/wbs3/working_dir/genomes/SDTS/SDTS.h1.chr.fasta
ref2=/media/nfs2/wbs/wbs3/working_dir/genomes/SDTS/SDTS.h2.chr.fasta

# scaffolding
ragtag.py scaffold -t 48 -o primary_ref1 ${ref1} $1.fasta
ragtag.py scaffold -t 48 -o primary_ref2 ${ref2} $1.fasta

