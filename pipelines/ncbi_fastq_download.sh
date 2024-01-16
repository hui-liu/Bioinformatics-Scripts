for i in $(cat list.txt)
do
echo $i
a=$(echo $i | cut -b 1-6)
b=$(echo $i | cut -b 10- | xargs printf "%03d")
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${a}/${b}/${i}/${i}_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${a}/${b}/${i}/${i}_2.fastq.gz
done
