samtools faidx ../$1.asm.bp.p_ctg.fa
awk '$2<5000000' ../$1.asm.bp.p_ctg.fa.fai | \
cut -f 1 | samtools faidx ../$1.asm.bp.p_ctg.fa -r - > Contig.fasta
blastn -query Contig.fasta \
-db /media/nfs2/wbs/wbs1/working_dir/genomes/Quercus_variabilis/CLM/organelle/organelle.fasta \
-out Contig_organelle.out -evalue 1e-5  -outfmt 6 -num_threads 24

python ~/bin/scripts/organelle_contigs.py \
Contig.fasta Contig_organelle.out > organelle.ids

python ~/bin/scripts/remove_seqs_by_ids.py ../$1.asm.bp.p_ctg.fa \
organelle.ids $1.asm.fa
