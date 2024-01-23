########
# OrthoFinder
#########
cd /media/nfs2/wbs/wbs1/working_dir/pan_genes/qva
orthofinder -f PEPs -t 48 -a 24 -X -o OrthoFinder

#
for i in $(ls PEPs/ | sed 's/.fasta//' | xargs)
do
grep ">" PEPs/${i}.fasta | sed "s/>/${i}\t/"
done > species_gene_map.txt

python ~/bin/scripts/qv_pangene/pangene_clusters_stats.py \
species_gene_map.txt \
OrthoFinder/Results_Jan16/Orthogroups/Orthogroups.txt 22

#########
# pangene
#########
cd /media/nfs2/wbs/wbs1/working_dir/pan_genes/qva/pan_gene
echo "SDTS CLM AH CB DJY DL FX GXNN HLT JJ LYG NX PG SJZ SS WM XA XR XY YT ZG ZS" | \
sed 's/ /\n/g' > order.txt

python ~/bin/scripts/qv_pangene/OrthoFinder2pangene.py \
../species_gene_map.txt \
order.txt \
../OrthoFinder/Results_Jan16/Orthogroups/Orthogroups.txt \
> pan_gene.csv


cut -d "," -f 3- pan_gene.csv | sed 's/,/\t/' | \
sed '1,1d' | sed 's/;/,/g' | \
awk '{gsub(",","\n"$1"\t");print $0}' | \
awk '{print $2"\t"$1}' | awk '$1!="NA"' > pan_gene_all_class.tsv


#########
# bootstrap
#########
cd /media/nfs2/wbs/wbs1/working_dir/pan_genes/qva/bootstrap
python ~/bin/scripts/qv_pangene/generating_boostrap_base.py \
../species_gene_map.txt \
../pan_gene/order.txt \
../OrthoFinder/Results_Jan16/Orthogroups/Orthogroups.txt

for k in {1..10}
do
for j in {1..100}
do
n=0
for i in `ls *pan_id.csv|shuf`
do
n=$((n+1))
cat $i >> temp.$j.$k
sort -u temp.$j.$k | wc -l >> result.pan.$j.$k
if [ $n == 1 ]; then
cat temp.$j.$k | wc -l >> result.core.$j.$k
fi
if [ $n == 2 ]; then
sort temp.$j.$k | uniq -c | awk '$1==2' | wc -l >> result.core.$j.$k
fi
if [ $n == 3 ]; then
sort temp.$j.$k | uniq -c | awk '$1==3' | wc -l >> result.core.$j.$k
fi
if [ $n == 4 ]; then
sort temp.$j.$k | uniq -c | awk '$1==4' | wc -l >> result.core.$j.$k
fi
if [ $n == 5 ]; then
sort temp.$j.$k | uniq -c | awk '$1==5' | wc -l >> result.core.$j.$k
fi
if [ $n == 6 ]; then
sort temp.$j.$k | uniq -c | awk '$1==6' | wc -l >> result.core.$j.$k
fi
if [ $n == 7 ]; then
sort temp.$j.$k | uniq -c | awk '$1==7' | wc -l >> result.core.$j.$k
fi
if [ $n == 8 ]; then
sort temp.$j.$k | uniq -c | awk '$1==8' | wc -l >> result.core.$j.$k
fi
if [ $n == 9 ]; then
sort temp.$j.$k | uniq -c | awk '$1==9' | wc -l >> result.core.$j.$k
fi
if [ $n == 10 ]; then
sort temp.$j.$k | uniq -c | awk '$1==10' | wc -l >> result.core.$j.$k
fi
if [ $n == 11 ]; then
sort temp.$j.$k | uniq -c | awk '$1==11' | wc -l >> result.core.$j.$k
fi
if [ $n == 12 ]; then
sort temp.$j.$k | uniq -c | awk '$1==12' | wc -l >> result.core.$j.$k
fi
if [ $n == 13 ]; then
sort temp.$j.$k | uniq -c | awk '$1==13' | wc -l >> result.core.$j.$k
fi
if [ $n == 14 ]; then
sort temp.$j.$k | uniq -c | awk '$1==14' | wc -l >> result.core.$j.$k
fi
if [ $n == 15 ]; then
sort temp.$j.$k | uniq -c | awk '$1==15' | wc -l >> result.core.$j.$k
fi
if [ $n == 16 ]; then
sort temp.$j.$k | uniq -c | awk '$1==16' | wc -l >> result.core.$j.$k
fi
if [ $n == 17 ]; then
sort temp.$j.$k | uniq -c | awk '$1==17' | wc -l >> result.core.$j.$k
fi
if [ $n == 18 ]; then
sort temp.$j.$k | uniq -c | awk '$1==18' | wc -l >> result.core.$j.$k
fi
if [ $n == 19 ]; then
sort temp.$j.$k | uniq -c | awk '$1==19' | wc -l >> result.core.$j.$k
fi
if [ $n == 20 ]; then
sort temp.$j.$k | uniq -c | awk '$1==20' | wc -l >> result.core.$j.$k
fi
if [ $n == 21 ]; then
sort temp.$j.$k | uniq -c | awk '$1==21' | wc -l >> result.core.$j.$k
fi
if [ $n == 22 ]; then
sort temp.$j.$k | uniq -c | awk '$1==22' | wc -l >> result.core.$j.$k
fi
done
rm temp.$j.$k
done &
done

for k in {1..10}
do
for j in {1..100}
do
cat result.core.$j.$k | xargs >> qva_core_genome_1000bs.txt
cat result.pan.$j.$k | xargs >> qva_pan_genome_1000bs.txt
done
done


#########
# kaks
#########
cd /media/nfs2/wbs/wbs1/working_dir/pan_genes/qva/kaks

ls ~/working_dir/gene_prediction/Quercus_variabilis/*/final/*PEP.longest.fasta | \
grep -vE "CLM.h2|SDTS.h2|SDTS.PEP" | xargs | \
awk '{print "cat",$0, "> PEP.fasta"}' | sh

ls ~/working_dir/gene_prediction/Quercus_variabilis/*/final/*CDS.longest.fasta | \
grep -vE "CLM.h2|SDTS.h2|SDTS.CDS" | xargs | \
awk '{print "cat",$0, "> CDS.fasta"}' | sh

#
awk -F "," '{print $3}' ../pan_gene/pan_gene.csv | sort | uniq -c
# 15909 core
# 18929 private
# 20989 shell
# 3411 softcore


group=../OrthoFinder/Results_Jan16/Orthogroups/Orthogroups.tsv
dos2unix ${group}
head -n 1 ${group} > Orthogroups.tsv

awk -F "," '$3=="core"{print $1}' ../pan_gene/pan_gene.csv | \
shuf -n 1591 | \
awk 'NR==FNR{a[$1]=$0;next}{$1 in a}{print a[$1]}' ${group} - \
>> Orthogroups.tsv

awk -F "," '$3=="softcore"{print $1}' ../pan_gene/pan_gene.csv | \
shuf -n 341 | \
awk 'NR==FNR{a[$1]=$0;next}{$1 in a}{print a[$1]}' ${group} - \
>> Orthogroups.tsv

awk -F "," '$3=="shell"{print $1}' ../pan_gene/pan_gene.csv | \
shuf -n 2099 | \
awk 'NR==FNR{a[$1]=$0;next}{$1 in a}{print a[$1]}' ${group} - \
>> Orthogroups.tsv

python ~/bin/scripts/pangenes/combinations4kaks_random.py \
Orthogroups.tsv > GenePairs.txt

awk -F "," '{print $1"\t"$3}' ../pan_gene/pan_gene.csv  | sed '1,1d' | \
awk 'NR==FNR{a[$1]=$1;next}{$1 in a}{print $0"\t"a[$1]}' Orthogroups.tsv - | \
awk 'NF==3' | cut -f 1,2 > cluster_class.tsv

#
cd /media/nfs2/wbs/wbs1/working_dir/pan_genes/qva/kaks
mkdir -p /dev/shm/KaKs_Calculator/OUT
python ~/bin/scripts/pangenes/generate_cal_kaks_list.py \
GenePairs.txt /dev/shm/KaKs_Calculator/OUT > chunks.txt

#
n=0
while read line
do
n=$[n+1]
echo "if [ \$SGE_TASK_ID -eq $n ];then"
echo "python ~/bin/scripts/pangenes/cal_kaks_YN.py $line PEP.fasta CDS.fasta /dev/shm"
echo "fi"
done < chunks.txt > cal_kaks.sh
qsub -tc 50 -t 1-$n -p -100 -l hostname=xnode03 -o cal_kaks.out -terse -cwd -j y -V -S /bin/bash cal_kaks.sh

find /dev/shm/KaKs_Calculator/OUT -name "*kaks" | xargs cat > GenePairs.kaks
sed -i '2,$ {/Sequence/d}' GenePairs.kaks

rm -r /dev/shm/KaKs_Calculator


#######
# pangene domian
#######
cd /mnt/Data_disk/liuhui/Qvariabilis/pan_gene/domain

ls ~/Qvariabilis/annotation/*/function/*interproscan.tsv | \
grep -vE "CLM.h2|SDTS.h2|Qa" | xargs | \
awk '{print "cat",$0, "| gzip -c > interproscan.tsv.gz"}' | sh

python ~/bin/scripts/pangene_domain.py \
interproscan.tsv.gz \
../pan_gene.csv \
> pangene_domain.tsv

#
zcat interproscan.tsv.gz | awk '$4=="Pfam"' | \
cut -f 1,5,6 > all_pfam.tsv
zcat interproscan.tsv.gz | cut -f 1,5,6 > all_interpro.tsv


#######
# GO
#######
cd /mnt/Data_disk/liuhui/Qvariabilis/pan_gene/GO
ls ~/Qvariabilis/annotation/*/function/*GO.tsv | \
grep -vE "CLM.h2|SDTS.h2|Qa" | xargs | \
awk '{print "cat",$0, "> all_GO.tsv"}' | sh

#######
# KEGG
#######
cd /mnt/Data_disk/liuhui/Qvariabilis/pan_gene/KEGG

ls ~/Qvariabilis/annotation/*/function/*KEGG.tsv | \
grep -vE "CLM.h2|SDTS.h2|Qa" | xargs | \
awk '{print "cat",$0, " | cut -f 1,2,4 > all_KEGG.tsv"}' | sh



#######
# expression
#######

cd /mnt/Data_disk/liuhui/Qvariabilis/expression
dir=/mnt/Data_disk/liuhui/data/Quercus_variabilis/RNAseq/X101SC23092165-Z01-J002/00.CleanData
for sample in AH CB DJY DL FX HLT JJ LYG NX PG SJZ WM XA XR XY YT ZG ZS
do
echo ${sample}
mkdir -p ${sample}
cd ${sample}
ln -s ${dir}/${sample}/*clean.fq.gz .
ls *gz | sed 's/_.*//' | uniq > sample_list
ln -s ~/Qvariabilis/genomes/${sample}/${sample}.fasta
ln -s ~/Qvariabilis/annotation/${sample}/${sample}.gene.longest.gtf
hisat2-build ${sample}.fasta ${sample}.fasta
sh ~/bin/scripts/call_TPM.sh \
${sample}.fasta \
${sample}.gene.longest.gtf \
sample_list
cd ..
done

# SDTS
dir=~/data/Quercus_variabilis/RNAseq/SDTS/
sample=SDTS
echo ${sample}
mkdir -p ${sample}
cd ${sample}
ln -s ${dir}/leaf_*.clean.fq.gz .
ls *gz | sed 's/_.*//' | uniq > sample_list
ln -s ~/Qvariabilis/genomes/${sample}/SDTS.h1.fasta
ln -s ~/Qvariabilis/annotation/${sample}/SDTS.h1.gene.longest.gtf
hisat2-build SDTS.h1.fasta SDTS.h1.fasta
sh ~/bin/scripts/call_TPM.sh \
SDTS.h1.fasta \
SDTS.h1.gene.longest.gtf \
sample_list


#
cd /mnt/Data_disk/liuhui/Qvariabilis/expression
for sample in AH CB DJY DL FX HLT JJ LYG NX PG SJZ WM XA XR XY YT ZG ZS
do
echo ${sample}
cd ${sample}
awk '$3=="transcript"' ${sample}.gene.longest.gtf | \
awk '{print $10"\t"$12}' | sed 's/"//g;s/;//g' > id_map.txt
python ~/bin/scripts/merge_tpm2.py id_map.txt > ${sample}_TPM.tsv
cd ..
done

# SDTS
cd /mnt/Data_disk/liuhui/Qvariabilis/expression/SDTS
awk '$3=="transcript"' SDTS.h1.gene.longest.gtf | \
awk '{print $10"\t"$12}' | sed 's/"//g;s/;//g' > id_map.txt
python ~/bin/scripts/merge_tpm2.py id_map.txt | sed 's/leaf/SDTS/' > SDTS_TPM.tsv


#
cd /mnt/Data_disk/liuhui/Qvariabilis/expression
cat */*TPM.tsv | sed '/gene_id/d' | sed '1igeneid\tTPM' > qv_leaf_TPM.tsv
