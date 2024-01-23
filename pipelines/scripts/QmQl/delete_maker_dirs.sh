cd /media/nfs2/wbs/wbs1/working_dir/gene_prediction/QmQl/$1/maker/$2
cat /media/nfs2/wbs/wbs1/working_dir/EDTA/QmQl/$1/$1.ids | \
awk '{print "rm -r", $1}' | sh

