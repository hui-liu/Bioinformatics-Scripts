source ~/.bash_conda
conda activate eggnog-mapper

emapper.py -i $1 \
-o $2 --target_taxa Viridiplantae --cpu 64

gzip $2.emapper.annotations
/media/nfs2/wbs/wbs1/bin/exe/python ~/bin/scripts/eggnog_mapper_representative_genes.py \
$3 $2.emapper.annotations.gz eggnog.tsv

cut -f 1,10 eggnog.tsv | awk '$2!="-"' | sed '1,1d' | \
awk '{gsub(/,/, "\n"$1"\t");print}' > eggnog_GO.tsv

