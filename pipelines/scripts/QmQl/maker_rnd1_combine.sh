# combine
awk '{print FILENAME"\t"$0}' */*/*master_datastore_index.log | \
sed 's#/#\t#2' | awk '{print $3"\t"$1"/"$4"\t"$5}' > master_datastore_index.log

fasta_merge -o rnd1 -d master_datastore_index.log
gff3_merge -n -s -d master_datastore_index.log > rnd1_noSeq.gff
# transcript alignments
awk '{ if ($2 == "est2genome") print $0 }' rnd1_noSeq.gff > rnd1.est2genome.gff
# protein alignments
awk '{ if ($2 == "protein2genome") print $0 }' rnd1_noSeq.gff > rnd1.protein2genome.gff

