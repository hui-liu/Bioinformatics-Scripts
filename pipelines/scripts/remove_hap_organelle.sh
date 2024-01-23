samtools faidx ../$1.asm.bp.$2.p_ctg.fa
awk '$2<5000000' ../$1.asm.bp.$2.p_ctg.fa.fai | \
cut -f 1 | samtools faidx ../$1.asm.bp.$2.p_ctg.fa -r - > $2.Contig.fasta
blastn -query $2.Contig.fasta \
-db /media/nfs2/wbs/wbs1/working_dir/genomes/Quercus_variabilis/CLM/organelle/organelle.fasta \
-out $2.Contig_organelle.out -evalue 1e-5  -outfmt 6 -num_threads 24

python ~/bin/scripts/organelle_contigs.py \
$2.Contig.fasta $2.Contig_organelle.out > $2.organelle.ids

python ~/bin/scripts/remove_seqs_by_ids.py ../$1.asm.bp.$2.p_ctg.fa \
$2.organelle.ids $1.$2.asm.fa
