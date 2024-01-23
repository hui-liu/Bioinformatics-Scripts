samtools faidx ../depth/$1.$2.asm.filtered.fa
awk '$2<5000000' ../depth/$1.$2.asm.filtered.fa.fai | \
cut -f 1 | samtools faidx ../depth/$1.$2.asm.filtered.fa -r - > $2.Contig.fasta
blastn -query $2.Contig.fasta \
-db /media/nfs2/wbs/wbs2/working_dir/organelle/Quercus_acutissima/Qa_organelle.fasta \
-out $2.Contig_organelle.out -evalue 1e-5  -outfmt 6 -num_threads 24

python ~/bin/scripts/organelle_contigs.py \
$2.Contig.fasta $2.Contig_organelle.out > $2.organelle.ids

python ~/bin/scripts/remove_seqs_by_ids.py ../depth/$1.$2.asm.filtered.fa \
$2.organelle.ids $1.$2.asm.fa
