source ~/.bash_conda
source activate EDTA_210

RepeatMasker -pa 64 -q -div 40 \
-cutoff 225 \
-gff \
-lib $2 \
$1.mod

perl -i -nle 's/\s+DNA\s+/\tDNA\/unknown\t/; print $_' $1.mod.out

perl /media/nfs2/wbs/wbs1/bin/anaconda3/envs/EDTA_210/bin/EDTA.pl \
--genome $1 -t 64 --step final --anno 1 \
--curatedlib $2 \
--rmout $1.mod.out
