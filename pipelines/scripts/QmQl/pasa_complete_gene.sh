genome=$1
gff3=$2
gff3_file_to_proteins.pl $gff3 $genome CDS  > evm.pasa.CDS.rnd2.raw.fasta

python ~/bin/scripts/get_complete_gene_gff.py evm.pasa.CDS.rnd2.raw.fasta \
$gff3 evm.pasa.rnd2.gff3

