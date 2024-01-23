n=0
while read line
do
n=$[n+1]
echo "if [ \$SGE_TASK_ID -eq $n ];then"
echo "$line"
echo "fi"
done < EVM_commands.list > run_EVM.sh
