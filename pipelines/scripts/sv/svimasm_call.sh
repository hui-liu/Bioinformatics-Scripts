# $1: sample
# $2: ref
source ~/.bash_conda
conda activate svimasm_env

svim-asm diploid --query_names --interspersed_duplications_as_insertions --min_sv_size 40 \
--sample $1 results $1.hap1.bam $1.hap2.bam $2
