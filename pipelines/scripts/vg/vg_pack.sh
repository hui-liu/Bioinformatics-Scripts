# $1: xg file
# $2: gam file
# $3: sample

vg pack \
-x $1 \
-g $2 \
-Q 5 \
-o $3.pack \
-t 48

