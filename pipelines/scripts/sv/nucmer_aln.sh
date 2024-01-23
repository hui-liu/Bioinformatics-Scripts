# $1: ref
# $2: qry
# $3: sample

nucmer --maxmatch -c 500 -b 500 -l 100 -t 48 $1 $2 --prefix $3
delta-filter -m -i 90 -l 100 $3.delta > $3.filtered.delta
