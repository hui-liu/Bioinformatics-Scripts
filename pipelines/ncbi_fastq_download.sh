for i in $(cat list.txt)
do
echo $i
a=$(echo $i | cut -b 1-6)
# b=(echo(echo i | cut -b 10- | xargs printf "%03d")
b=$(echo $i | cut -b 10-)
c=$(printf "%03d" "${b#0}")
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${a}/${b}/${i}/${i}_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${a}/${b}/${i}/${i}_2.fastq.gz
done
