samtools faidx ../depth/$1.asm.filtered.fa
awk '$2<5000000' ../depth/$1.asm.filtered.fa.fai | \
cut -f 1 | samtools faidx ../depth/$1.asm.filtered.fa -r - > Contig.fasta
blastn -query Contig.fasta \
-db /media/nfs2/wbs/wbs2/working_dir/organelle/Quercus_acutissima/Qa_organelle.fasta \
-out Contig_organelle.out -evalue 1e-5  -outfmt 6 -num_threads 24

python ~/bin/scripts/organelle_contigs.py \
Contig.fasta Contig_organelle.out > organelle.ids

python ~/bin/scripts/remove_seqs_by_ids.py ../depth/$1.asm.filtered.fa \
organelle.ids $1.asm.fa
