samtools faidx $1/ragtag.scaffold.fasta
grep Chr $1/ragtag.scaffold.fasta.fai | cut -f 1 | sort -V | \
samtools faidx $1/ragtag.scaffold.fasta -r - | \
sed 's/_RagTag//' | sed '/>/s/[ab]$//' > $2.$3.chr.fasta
