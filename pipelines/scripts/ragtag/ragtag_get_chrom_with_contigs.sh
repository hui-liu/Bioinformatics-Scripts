#samtools faidx $1/ragtag.scaffold.fasta
grep Chr $1/ragtag.scaffold.fasta.fai | cut -f 1 | sort -V > ids
grep -v Chr $1/ragtag.scaffold.fasta.fai | cut -f 1 | sort -V >> ids

samtools faidx $1/ragtag.scaffold.fasta -r ids | \
sed 's/_RagTag//' | sed '/>/s/[ab]$//' | \
sed "/Chr/s/>/>$2./" > $2.fasta
rm ids
