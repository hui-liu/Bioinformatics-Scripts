# $1: xg file
# $2: snarls file
# $3: sample

vg call \
$1 \
-r $2 \
-s $3 \
-k $3.pack \
-t 48 \
-a | bgzip -c > $3.vcf.gz

