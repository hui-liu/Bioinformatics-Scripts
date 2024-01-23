n=0
for i in $(cat $1)
do
n=$[n+1]
mkdir -p $i/tmp
cp ~/bin/scripts/QmQl/maker/rnd3/*.ctl $i
sed -i "s#^genome=#&$2/$i.fa#" $i/maker_opts.ctl
sed -i "s#^est_gff=#&$3#" $i/maker_opts.ctl
sed -i "s#^protein_gff=#&$4#" $i/maker_opts.ctl
sed -i "s#^rm_gff=#&$5#" $i/maker_opts.ctl
sed -i "s#^snaphmm=#&$6#" $i/maker_opts.ctl
sed -i "s#^gmhmm=#&$7#" $i/maker_opts.ctl
sed -i "s#^augustus_species=#&$8#" $i/maker_opts.ctl
sed -i "s#^TMP=#&$9/$i/tmp#" $i/maker_opts.ctl
echo "if [ \$SGE_TASK_ID -eq $n ];then"
echo "cd $9/$i"
echo "source ~/.bash_conda"
echo "source activate maker"
echo "export PATH=/media/nfs2/wbs/wbs1/bin/Augustus/bin:\$PATH"
echo "export PATH=/media/nfs2/wbs/wbs1/bin/Augustus/scripts:\$PATH"
echo "export AUGUSTUS_CONFIG_PATH=/media/nfs2/wbs/wbs1/bin/Augustus/config/"
echo "maker -cpus 16 --ignore_nfs_tmp -base $i"
echo "fi"
done > $8_rnd3.sh
