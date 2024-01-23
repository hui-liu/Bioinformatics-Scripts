source ~/.bash_conda
conda activate ragtag

ref=/media/nfs2/wbs/wbs1/working_dir/genomes/Quercus_mongolica2_update/QuMo.fa

# scaffolding
ragtag.py scaffold -t 48 -o result ${ref} $1.fasta

