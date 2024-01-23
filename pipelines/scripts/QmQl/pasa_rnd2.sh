sed "s#../../evm/evm.all.gff3#../rnd1/pasa.sqlite.gene_structures_post_PASA_updates.$1.gff3#" \
../rnd1/run_annotationCompare_rnd1.sh | head -n 1 > run_annotationCompare_rnd2.sh

