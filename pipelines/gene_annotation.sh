#########
# EDTA
#########
genomedir=/media/nfs2/wbs/wbs1/working_dir/genomes/QmQl
EDTAdir=/media/nfs2/wbs/wbs1/working_dir/EDTA/QmQl
cd ${EDTAdir}

n=0
for i in $(ls ~/working_dir/genomes/QmQl/ | xargs)
do
n=$[n+1]
mkdir -p ${i}
echo "if [ \$SGE_TASK_ID -eq $n ];then"
echo "cd ${EDTAdir}/${i}"
echo "/media/nfs2/wbs/wbs1/bin/scripts/annotation/run_EDTA.sh \
${genomedir}/${i}/${i}.fasta &> ${i}_EDTA.log"
echo "fi"
done > run_EDTA.sh
qsub -pe mpi 1 -tc 2 -t 1-$n -o run_EDTA.out -terse -cwd -j y -V -S /bin/bash run_EDTA.sh

#
EDTAdir=/media/nfs2/wbs/wbs1/working_dir/EDTA/QmQl
cd ${EDTAdir}

for i in $(ls ~/working_dir/genomes/QmQl/ | xargs)
do
cd ${i}
python2 ~/bin/scripts/edta_repeat2gff3.py \
${i}.fasta.mod.EDTA.anno/${i}.fasta.mod.out ${i}.maker.gff3
cd ..
done

########
# split genome 
######## 
#split the genome by chromosome and all of the contigs were put into Contigs.fa
EDTAdir=/media/nfs2/wbs/wbs1/working_dir/EDTA/QmQl
cd ${EDTAdir}

for i in $(ls ~/working_dir/genomes/QmQl/ | xargs)
do
cd ${i}
python2 ~/bin/scripts/split_contigs.py \
${i}.fasta.mod.MAKER.masked 1500000 > ${i}.ids
cd ..
done

###########
# make the repeatmasker gff3 file for maker
###########
dir=/media/nfs2/wbs/wbs1/working_dir/repeatmasker/QmQl
EDTAdir=/media/nfs2/wbs/wbs1/working_dir/EDTA/QmQl
cd ${dir}

n=0
for i in $(ls ~/working_dir/genomes/QmQl/ | xargs)
do
n=$[n+1]
mkdir -p ${i}/result
echo ">(N)n#Dummy_repeat @root  [S:25]
nnnnnnnnnnnnnnnnnn" > ${i}/simple.lib
echo "if [ \$SGE_TASK_ID -eq $n ];then"
echo "cd ${dir}/${i}"
echo "/media/nfs2/wbs/wbs1/bin/RepeatMasker-open-4-0-7/RepeatMasker \
${EDTAdir}/${i}/${i}.fasta.mod.MAKER.masked \
-dir result \
-pa 48 \
-lib simple.lib \
-gff"
echo "/media/nfs2/wbs/wbs1/bin/exe/python ~/bin/scripts/repeat2gff3.py \
result/${i}.fasta.mod.MAKER.masked.out result/${i}.fasta.mod.MAKER.masked.out.gff3"
echo "cat ../../../EDTA/QmQl/${i}/${i}.maker.gff3 \
result/${i}.fasta.mod.MAKER.masked.out.gff3 > ${i}.repeat.gff3"
echo "fi"
done > run_repeatmasker.sh

qsub -pe mpi 1 -tc 2 -t 1-$n -o run_repeatmasker.out -terse -cwd -j y -V -S /bin/bash run_repeatmasker.sh


#########
# BUSCO
#########
genomedir=/media/nfs2/wbs/wbs1/working_dir/genomes/QmQl
BUSCOdir=/media/nfs2/wbs/wbs1/working_dir/BUSCO/QmQl
cd ${BUSCOdir}

n=0
for i in $(ls ~/working_dir/genomes/QmQl/ | xargs)
do
n=$[n+1]
mkdir -p ${i}
echo "if [ \$SGE_TASK_ID -eq $n ];then"
echo "cd ${BUSCOdir}/${i}"
echo "/media/nfs2/wbs/wbs1/bin/scripts/annotation/run_busco.sh \
${genomedir}/${i}/${i}.fasta ${i} geno &> ${i}_busco.log"
echo "fi"
done > run_busco.sh
qsub -pe mpi 1 -tc 2 -t 1-$n -o run_busco.out -terse -cwd -j y -V -S /bin/bash run_busco.sh

#########
# AUGUSTUS
#########
cd /media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl
n=0
for sample in $(ls ~/working_dir/genomes/QmQl/ | xargs)
do
n=$[n+1]
dir=/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/AUGUSTUS
BUSCO=/media/nfs2/wbs/wbs1/working_dir/BUSCO/QmQl/${sample}/${sample}
genome=/media/nfs2/wbs/wbs1/working_dir/genomes/QmQl/${sample}/${sample}.fasta
echo "if [ \$SGE_TASK_ID -eq $n ];then"
echo "mkdir -p ${dir}"
echo "cd ${dir}"
echo "~/bin/scripts/AUGUSTUS_train.sh ${sample} ${BUSCO} ${genome}"
echo "fi"
done > run_AUGUSTUS.sh
qsub -pe mpi 1 -tc 1 -t 1-$n -o run_AUGUSTUS.out -terse -cwd -j y -V -S /bin/bash run_AUGUSTUS.sh

#########
# GeneMark
#########
cd /media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl
n=0
for sample in $(ls ~/working_dir/genomes/QmQl/ | xargs)
do
n=$[n+1]
dir=/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/GeneMark
genome=/media/nfs2/wbs/wbs1/working_dir/EDTA/QmQl/${sample}/${sample}.fasta.mod.MAKER.masked
echo "if [ \$SGE_TASK_ID -eq $n ];then"
echo "mkdir -p ${dir}"
echo "cd ${dir}"
echo "~/bin/scripts/GeneMark.sh ${genome}"
echo "fi"
done > run_GeneMark.sh

qsub -l hostname=xnode01 -pe mpi 1 -tc 3 -t 1-$n -o run_GeneMark.out -terse -cwd -j y -V -S /bin/bash run_GeneMark.sh


########
# maker
########

cd /media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl
# rnd1
for sample in $(ls ~/working_dir/genomes/QmQl/ | xargs)
do
dir=/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/maker
mkdir -p ${dir}
cd ${dir}
mkdir rnd1
mkdir rnd2
mkdir rnd3
cd ${dir}/rnd1
~/bin/scripts/QmQl/maker_rnd1.sh \
/media/nfs2/wbs/wbs1/working_dir/EDTA/QmQl/${sample}/${sample}.ids \
/media/nfs2/wbs/wbs1/working_dir/EDTA/QmQl/${sample} \
/media/nfs2/wbs/wbs1/working_dir/rnaseq/QmQl/${sample}/${sample}_transcriptome.fasta \
/media/nfs2/wbs/wbs1/working_dir/repeatmasker/QmQl/${sample}/${sample}.repeat.gff3 \
${dir}/rnd1 \
${sample}
done

# rnd2
for sample in $(ls ~/working_dir/genomes/QmQl/ | xargs)
do
dir=/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/maker
cd ${dir}/rnd2
~/bin/scripts/QmQl/maker_rnd2.sh \
/media/nfs2/wbs/wbs1/working_dir/EDTA/QmQl/${sample}/${sample}.ids \
/media/nfs2/wbs/wbs1/working_dir/EDTA/QmQl/${sample} \
/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/maker/rnd1/rnd1.est2genome.gff \
/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/maker/rnd1/rnd1.protein2genome.gff \
/media/nfs2/wbs/wbs1/working_dir/repeatmasker/QmQl/${sample}/${sample}.repeat.gff3 \
/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/snap/rnd1/genome.hmm \
/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/GeneMark/gmhmm.mod \
${sample} \
${dir}/rnd2
done

# rnd3
for sample in $(ls ~/working_dir/genomes/QmQl/ | xargs)
do
dir=/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/maker
cd ${dir}/rnd3
~/bin/scripts/QmQl/maker_rnd3.sh \
/media/nfs2/wbs/wbs1/working_dir/EDTA/QmQl/${sample}/${sample}.ids \
/media/nfs2/wbs/wbs1/working_dir/EDTA/QmQl/${sample} \
/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/maker/rnd2/rnd2.est2genome.gff \
/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/maker/rnd2/rnd2.protein2genome.gff \
/media/nfs2/wbs/wbs1/working_dir/repeatmasker/QmQl/sample/{sample}/{sample}.repeat.gff3 \
/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/snap/rnd2/genome.hmm \
/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/GeneMark/gmhmm.mod \
${sample} \
${dir}/rnd3
done

######
# An example
#######
sample=DD1Qm
cd /media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/maker/rnd1
n=(grepSGETASKID(grep SGE_TASK_ID  {sample}_rnd1.sh | wc -l)
qsub -pe mpi 1 -tc 25 -t 1-n−on -o {sample}_rnd1.out -terse -cwd -j y -V -S /bin/bash ${sample}_rnd1.sh

~/bin/scripts/QmQl/maker_rnd1_combine.sh

dir=/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/snap/rnd1
mkdir -p ${dir}
cd ${dir}
~/bin/scripts/QmQl/snap_rnd1.sh

sample=DD1Qm
~/bin/scripts/QmQl/delete_maker_dirs.sh ${sample} rnd1

sample=DD1Qm
cd /media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/maker/rnd2
n=$(grep SGE_TASK_ID  ${sample}_rnd2.sh | wc -l)
qsub -pe mpi 1 -tc 25 -t 1-$n -o ${sample}_rnd2.out -terse -cwd -j y -V -S /bin/bash ${sample}_rnd2.sh

~/bin/scripts/QmQl/maker_rnd2_combine.sh

dir=/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/snap/rnd2
mkdir -p ${dir}
cd ${dir}
~/bin/scripts/QmQl/snap_rnd2.sh

sample=DD1Qm
~/bin/scripts/QmQl/delete_maker_dirs.sh ${sample} rnd2

sample=DD1Qm
cd /media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/maker/rnd3
n=$(grep SGE_TASK_ID  ${sample}_rnd3.sh | wc -l)
qsub -pe mpi 1 -tc 25 -t 1-$n -o ${sample}_rnd3.out -terse -cwd -j y -V -S /bin/bash ${sample}_rnd3.sh

~/bin/scripts/QmQl/maker_rnd3_combine.sh

sample=DD1Qm
~/bin/scripts/QmQl/delete_maker_dirs.sh ${sample} rnd3

sample=DD1Qm
dir=/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}
mkdir -p ${dir}/evm
cd ${dir}/evm
~/bin/scripts/QmQl/partition_EVM_inputs.sh \
${dir} \
/media/nfs2/wbs/wbs1/working_dir/genomes/QmQl/${sample}/${sample}.fasta

qqsub run_partition_EVM_inputs.sh

~/bin/scripts/QmQl/EVM.sh

n=$(cat EVM_commands.list| wc -l)
qsub -tc 128 -t 1-$n -o run_EVM.out -terse -cwd -j y -V -S /bin/bash run_EVM.sh

sample=DD1Qm
~/bin/scripts/QmQl/EVM_combine.sh \
/media/nfs2/wbs/wbs1/working_dir/genomes/QmQl/${sample}/${sample}.fasta

sample=DD1Qm
dir=/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}
mkdir -p ${dir}/pasa
cd ${dir}/pasa
ln ~/working_dir/rnaseq/QmQl/${sample}/${sample}_transcriptome.fasta

echo "seqclean ${sample}_transcriptome.fasta -v $HOME/bin/PASApipeline-v2.5.2/UniVec/UniVec" > run_seqclean.sh
qqsub run_seqclean.sh

sample=DD1Qm
dir=/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/pasa
cd ${dir}
mkdir rnd1 rnd2
~/bin/scripts/QmQl/pasa_align.sh \
${dir} \
/media/nfs2/wbs/wbs1/working_dir/genomes/QmQl/${sample}/${sample}.fasta \
${sample}

qqsub -l hostname=xnode01 run_alignAssembly.sh

~/bin/PASApipeline-v2.5.2/misc_utilities/pasa_gff3_validator.pl ../evm/evm.all.gff3

sample=DD1Qm
dir=/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/pasa
cd ${dir}/rnd1
~/bin/scripts/QmQl/pasa_rnd1.sh ${sample}
qqsub -l hostname=xnode01 run_annotationCompare_rnd1.sh

sample=DD1Qm
dir=/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/pasa
cd ${dir}/rnd2
num=$(ls ../rnd1/*gff3 | awk -F "." '{print $(NF-1)}')
~/bin/scripts/QmQl/pasa_rnd2.sh $num
qqsub -l hostname=xnode01 run_annotationCompare_rnd2.sh

num=$(ls *gff3 | awk -F "." '{print $(NF-1)}')
~/bin/scripts/QmQl/pasa_complete_gene.sh \
/media/nfs2/wbs/wbs1/working_dir/genomes/QmQl/${sample}/${sample}.fasta \
pasa.sqlite.gene_structures_post_PASA_updates.$num.gff3

sample=DD1Qm
dir=/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/filter
mkdir -p ${dir}
cd ${dir}
python2 ~/bin/scripts/ann_redup.py \
../pasa/rnd2/evm.pasa.rnd2.gff3 \
evm.pasa.filtered.gff3 > filter.log

sample=DD1Qm
dir=/media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/final
mkdir -p ${dir}
cd ${dir}
~/bin/scripts/QmQl/rename.sh \
/media/nfs2/wbs/wbs1/working_dir/genomes/QmQl/${sample}/${sample}.fasta \
${sample}
