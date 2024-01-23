n=0
for i in $(cat $1)
do
n=$[n+1]
mkdir -p $i/tmp
cp ~/bin/scripts/QmQl/maker/rnd1/*.ctl $i
sed -i "s#^genome=#&$2/$i.fa#" $i/maker_opts.ctl
sed -i "s#^est=#&$3#" $i/maker_opts.ctl
sed -i "s#^rm_gff=#&$4#" $i/maker_opts.ctl
sed -i "s#^TMP=#&$5/$i/tmp#" $i/maker_opts.ctl
echo "if [ \$SGE_TASK_ID -eq $n ];then"
echo "cd $5/$i"
echo "source ~/.bash_conda"
echo "source activate maker"
echo "maker -cpus 16 --ignore_nfs_tmp -base $i"
echo "fi"
done > $6_rnd1.sh
