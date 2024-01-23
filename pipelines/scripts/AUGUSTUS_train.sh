species=$1
dir=$2
genome=$3
python=/media/nfs2/wbs/wbs1/bin/anaconda3/envs/py2/bin/python2.7
#
cat $dir/run_embryophyta_odb10/busco_sequences/single_copy_busco_sequences/*faa > single_copy_busco.faa

# removes genes with distances <=2000
# remove genes with identity >=0.8 on aa level
grep ">" single_copy_busco.faa | sed 's/>//' | awk '{print $1"\t"$1}' | \
sed 's/:/\t/;s/-/\t/' | sort -V > single_copy_busco.bed

bedtools closest -nonamecheck \
-a single_copy_busco.bed -b single_copy_busco.bed -k 2000 -N -d > single_copy_busco.dist
cd-hit -o single_copy_busco.faa.cdhit -c 0.8 -i single_copy_busco.faa -p 1 -d 0 -b 3 -T 8 -M 16000

cat <(grep -B 1 "^[1-9]" single_copy_busco.faa.cdhit.clstr | awk '{print $3}' | sed 's/>//' | sed 's/\.\.\.$//' | sed '/^$/d') \
<(awk '$9<=2000 && $9>=0' single_copy_busco.dist | cut -f 4,8 | sed 's/\t/\n/' | awk '!a[$1]++') | \
awk '!a[$1]++' > badgenes.txt

${python} ~/bin/scripts/remove_seqs_by_ids.py single_copy_busco.faa badgenes.txt single_copy_busco_good.faa

${python} ~/bin/scripts/get_busco_cds_gff.py single_copy_busco_good.faa \
$dir/run_embryophyta_odb10/full_table.tsv \
$dir/run_embryophyta_odb10/augustus_output/gff single_copy_busco_good.gff


# Convert GFF file to Genbank format file
gff2gbSmallDNA.pl single_copy_busco_good.gff $genome 1000 genes.gb

# Split gene structure set into training and test set
num=$(grep "LOCUS" genes.gb | wc -l)
num=$(echo $num -1000 | bc)
randomSplit.pl genes.gb $num

# CREATE A META PARAMETERS FILE FOR YOUR SPECIES
new_species.pl --species=$species

# Initial training
etraining --species=$species genes.gb.train

# predict the genes in genes.gb.train ab initio
augustus --species=$species genes.gb.test | tee firsttest.out
grep -A 22 Evaluation firsttest.out > firsttest.log

# RUN THE SCRIPT optimize_augustus.pl
optimize_augustus.pl --rounds=5 --cpus=64 --species=$species genes.gb.train
etraining --species=$species genes.gb.train
augustus --species=$species genes.gb.test | tee secondtest.out

grep -A 22 Evaluation secondtest.out > secondtest.log


