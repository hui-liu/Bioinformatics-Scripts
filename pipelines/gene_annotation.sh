sample=DD1Qm
cd /media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/${sample}/maker/rnd1
n=$(grep SGE_TASK_ID  ${sample}_rnd1.sh | wc -l)
qsub -pe mpi 1 -tc 25 -t 1-$n -o ${sample}_rnd1.out -terse -cwd -j y -V -S /bin/bash ${sample}_rnd1.sh

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
