genome=$1
sample=$2
samtools faidx $genome
cat $1.fai | \
awk '{print $1,$1,$1}' | sed "s/${sample}//3" | \
awk -v var=${sample} \
'{printf "%s\t%s\t"var"%03d\n",$1,$2,$3}' > id_conversion.tsv

~/bin/scripts/rename_gff.sh \
../filter/evm.pasa.filtered.gff3 \
id_conversion.tsv \
100 \
${sample}


python2 ~/bin/scripts/gff3_to_gtf.py ${sample}.gene.gff3 > ${sample}.gene.gtf

gff3_file_to_proteins.pl ${sample}.gene.gff3 $genome prot  > ${sample}.PEP.fasta
gff3_file_to_proteins.pl ${sample}.gene.gff3 $genome CDS  > ${sample}.CDS.fasta

# longest
python ~/bin/scripts/get_longest_alternative_seq.py ${sample}.gene.gff3 \
${sample}.CDS.fasta ${sample}.PEP.fasta \
${sample}.CDS.longest.fasta ${sample}.PEP.longest.fasta

python ~/bin/scripts/get_longest_feature_gff.py ${sample}.PEP.longest.fasta \
${sample}.gene.gff3  ${sample}.gene.longest.gff3

python2 ~/bin/scripts/gff3_to_gtf.py ${sample}.gene.longest.gff3 > ${sample}.gene.longest.gtf

#
sed '1i##gff-version\t3' ${sample}.gene.longest.gff3 | \
gt gff3 -tidy -addids no -retainids yes -addintrons yes | \
grep -v "^#" | python ~/bin/scripts/intron_attris.py ${sample}.gene.intron.longest.gff3
