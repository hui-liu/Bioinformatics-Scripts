cd /mnt/crick/data/swasp

#############################
# SET UP PATHS
#############################
rawgff='/mnt/crick/PLAZA/original_download/populus_tremula/v2.2/gff/Potra02_genes.gff'
gff='/mnt/crick/data/swasp/annotation/Potra02_genes.gff3'
fasta='/mnt/crick/PLAZA/original_download/populus_tremula/v2.2/fasta/Potra02_genome.fasta'
mygenome='/mnt/crick/data/swasp/annotation/Potra02.genome'
CDSbedout='/mnt/crick/data/swasp/annotation/Potra_CDS.bed'
fastaCDSout='/mnt/crick/data/swasp/annotation/Potra_CDS.tab'
degeneratebedout='/mnt/crick/data/swasp/annotation/Potra'
utr3='/mnt/crick/data/swasp/annotation/Potra_utr3.bed'
utr5='/mnt/crick/data/swasp/annotation/Potra_utr5.bed'
intron='/mnt/crick/data/swasp/annotation/Potra_intron.bed'
upstream='/mnt/crick/data/swasp/annotation/Potra_upstream.bed'
downstream='/mnt/crick/data/swasp/annotation/Potra_downstream.bed'
intergenic='/mnt/crick/data/swasp/annotation/Potra_intergenic.bed'

###############################
# add introns
###############################
gt gff3 -tidy -addids no -retainids yes -addintrons yes ${rawgff} | grep -v "^#" > ${gff}

##############################
# bed files: cds, utr3, utr5, intron, upstream, downstream, and intergenic sites
##############################

# get bed file of cds and shift to match phase
bash ~/bin/Degeneracy/gff2bed.sh ${gff} CDS | awk -f ~/bin/Degeneracy/gffphaseshift.awk  > ${CDSbedout}

# get bed file of three prime UTR
awk '$3=="three_prime_UTR"' ${gff} | gff2bed | awk -v OFS="\t" '{print $1,$2,$3,"utr3",$5,$6}' > ${utr3}

# get bed file of five prime UTR
awk '$3=="five_prime_UTR"' ${gff} | gff2bed | awk -v OFS="\t" '{print $1,$2,$3,"utr5",$5,$6}' > ${utr5}

# get bed file of intron
awk '$3=="intron"' ${gff} | gff2bed | awk -v OFS="\t" '{print $1,$2,$3,"intron",$5,$6}' > ${intron}

# get bed file of upstream
awk '$3=="gene"' ${gff} | bedtools flank -i - -g ${mygenome} -l 2000 -r 0 -s | \
awk -v OFS="\t" '{print $1,$4,$5,"upstream",$7,$8}' > ${upstream}

# get bed file of downstream
awk '$3=="gene"' ${gff} | bedtools flank -i - -g ${mygenome} -l 0 -r 2000 -s | \
awk -v OFS="\t" '{print $1,$4,$5,"downstream",$7,$8}' > ${downstream}

# get bed file of intergenic
awk '$3=="gene"' ${gff} | bedtools flank -i - -g ${mygenome} -b 2000 | \
awk '$9==prev {rec=rec "\t"$5}$9!=prev {if (NR>1) print rec; rec=$0}{prev=$9}END{print rec}' | \
awk -v OFS="\t" '$5=$10' | cut -f 1-9 | sort -V | \
bedtools complement -i - -g ${mygenome}| \
awk -v OFS="\t" '{print $1,$2,$3,"intergenic"}' > ${intergenic}

##############################
# 0 fold and 4 fold degenerate sites
##############################
# USE BED FILE AND FASTA TO GET FILE OF POS AND SEQUENCE
bedtools getfasta -s -tab -name -fi ${fasta} -bed ${CDSbedout} > ${fastaCDSout}
# CONVERT FASTA DNA SEQUENCE INTO CODONS, FLIP FOR PHASE, AND REPORT 4FOLD SITES
python ~/bin/Degeneracy/degeneracy.py -i ${fastaCDSout} -o ${degeneratebedout}

