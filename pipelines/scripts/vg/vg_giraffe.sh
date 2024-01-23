# $1: gbz
# $2: min
# $3: dist
# $4: path to fastq
# $5: sample

vg giraffe \
-Z $1 \
-m $2 \
-d $3 \
-x $4 \
-f $5/$6_filter_R1.fastq.gz \
-f $5/$6_filter_R2.fastq.gz \
-t 48 \
-p \
--read-group "ID:$6 SM:$6" \
--sample $6 \
> $6.gam

