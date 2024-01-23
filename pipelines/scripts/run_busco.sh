source ~/.bash_conda
conda activate busco
export PATH=$HOME/bin/Augustus/bin:$PATH
export PATH=$HOME/bin/Augustus/scripts:$PATH
export AUGUSTUS_CONFIG_PATH=$HOME/bin/Augustus/config/

# $1: input
# $2: output
# $3: geno or genome
#     tran or transcriptome
#     prot or proteins

/media/nfs2/wbs/wbs1/bin/anaconda3/envs/busco/bin/busco \
-i $1 \
-o $2 -m $3 --offline \
-l /media/nfs2/wbs/wbs1/bin/anaconda3/envs/busco/databases/embryophyta_odb10 \
-c 24 --augustus
