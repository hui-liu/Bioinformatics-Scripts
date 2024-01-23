# $1: input pep
source ~/.bash_conda
exec_annotation -o KofamKOALA.tsv --cpu 48 --format detail-tsv $1

gzip KofamKOALA.tsv
/media/nfs2/wbs/wbs1/bin/exe/python ~/bin/scripts/top_x_KofamKOALA.py KofamKOALA.tsv.gz 5 KofamKOALA_best5.tsv

gzip KofamKOALA_best5.tsv
/media/nfs2/wbs/wbs1/bin/exe/python ~/bin/scripts/KofamKOALA_representative_genes.py \
$2 \
KofamKOALA_best5.tsv.gz KEGG.tsv

