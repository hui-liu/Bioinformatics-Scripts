source ~/.bash_conda
source activate EDTA_210
perl /media/nfs2/wbs/wbs1/bin/anaconda3/envs/EDTA_210/bin/EDTA.pl \
--genome $1 \
--overwrite 1 --sensitive 1 --anno 1 --threads 64
