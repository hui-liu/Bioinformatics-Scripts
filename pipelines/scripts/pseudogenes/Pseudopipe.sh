# Prepare
pre=$1
genome=$2
gff=$3
pep=$4
dir=$5

#pre=QGR
#genome=/mnt/Data_disk/liuhui/genomes/Quercus_griffithi/QuGr.fa
#gff=/mnt/Data_disk/liuhui/genomes/Quercus_griffithi/Quercus_griffithi.gene.longest.gff3
#pep=/mnt/Data_disk/liuhui/genomes/Quercus_griffithi/Quercus_griffithi.PEP.longest.fasta
#dir=/mnt/Data_disk/liuhui/Quercus_pangenome/Pseudogene

# Input files configuration
cd $dir
mkdir -p {ppipe_output/$pre,ppipe_input/$pre}
mkdir -p ppipe_input/$pre/{dna,pep,mysql}
cp $genome ppipe_input/$pre/dna/
cp $pep ppipe_input/$pre/pep/

cd $dir/ppipe_input/$pre/dna/
# a list file for all unmasked dna divided into different chromosomes in FASTA format
python /media/nfs2/wbs/wbs1/bin/scripts/pseudogenes/split_seqs_by_chr.py $genome $pre
cd $dir/ppipe_input/$pre/mysql/
python /media/nfs2/wbs/wbs1/bin/scripts/pseudogenes/extract_exon_by_chr.py $gff


# Pseudopipe configuration
PSEUDOPIPE_HOME=/media/nfs2/wbs/wbs1/bin/pgenes/pseudopipe
export pseudopipe=$PSEUDOPIPE_HOME/core/runScripts.py
export genPgeneResult=$PSEUDOPIPE_HOME/ext/genPgeneResult.sh
export genFullAln=$PSEUDOPIPE_HOME/ext/genFullAln.sh
export fastaSplitter=$PSEUDOPIPE_HOME/ext/splitFasta.py
export sqDummy=$PSEUDOPIPE_HOME/ext/sqDummy.py
export blastHandler=$PSEUDOPIPE_HOME/core/processBlastOutput.py
export extractExLoc=$PSEUDOPIPE_HOME/core/extractKPExonLocations.py

# Alignment tools configuration
export formatDB=/media/nfs2/wbs/wbs1/bin/blast-2.2.25/bin/formatdb
export blastExec=/media/nfs2/wbs/wbs1/bin/blast-2.2.25/bin/blastall
export fastaExec=/media/nfs2/wbs/wbs1/bin/fasta34/tfasty34

# Python configuration
export pythonExec=python


outDir=$dir/ppipe_output/$pre
rmkDir=$dir/ppipe_input/$pre/dna/${genome##*/}
dnaTmp=$dir/ppipe_input/$pre/dna/$pre.%s.fa
pepDir=$dir/ppipe_input/$pre/pep/${pep##*/}
emkTmp=$dir/ppipe_input/$pre/mysql/%s_exLocs

#
inputDNA=$outDir/dna/dna_all.fa
inputPEP=$outDir/pep/pep_all.fa

#
jobsExec=$sqDummy

# (1) Making directories
cd $outDir
mkdir dna pep
#
mkdir -p blast/{stamps,output,status,processed}

#
mkdir -p pgenes/minus/{log,stamp}
mkdir -p pgenes/plus/{log,stamp}

# (2) Copying sequences
cp $rmkDir $inputDNA
cp $pepDir $inputPEP

# (3) Fomatting the DNAs
cd dna
$formatDB -i $inputDNA -o T -p F

# (4) Preparing the blast jobs
cd ../blast

pepNum=`grep -c '>' $inputPEP`
$pythonExec $fastaSplitter $inputPEP $((($pepNum+399)/400)) 'split%04d'
for f in `ls split*` 
do
	echo "( cd $(pwd); touch stamps/${f}.Start ; ( $blastExec -p tblastn -m 8 -z 4.7e8 -e .1 -d $inputDNA -i $f -o output/${f}.Out ; touch stamps/${f}.Stamp ) >status/${f}.Status 2>&1 )" 
done > jobs

#$pythonExec $jobsExec jobs
ParaFly -c jobs -CPU 96


# (5) Processing blast output
cd processed

echo "\"(cd $(pwd); $pythonExec $blastHandler $inputPEP  'split\d{4}.Out\Z' ../output; touch processed.stamp)\"" > jobs
$pythonExec $jobsExec jobs

# (6) Running Pseudopipe on both strands
cd ../../pgenes

for t in 'M' 'P'
do
	echo 'Working on '$t' strand'

	if [ $t = 'M' ]
	then
		cd minus
	else
		cd plus
	fi

	echo "export BlastoutSortedTemplate=${outDir}/blast/processed/%s_${t}_blastHits.sorted;export ChromosomeFastaTemplate=${dnaTmp};export ExonMaskTemplate=${emkTmp};export ExonMaskFields='2 3';export FastaProgram=${fastaExec};export ProteinQueryFile=${inputPEP}" > setenvPipelineVars
    source setenvPipelineVars
	ms=$(cd ../../blast/processed ; for f in *_${t}_*sorted; do echo ${f/_${t}_blastHits.sorted/}; done)
    # enumerate each chromosome
	for c in $ms
	do
        	echo "\"(cd $(pwd); touch stamp/$c.Start ; $pythonExec $pseudopipe $c > log/$c.log 2>&1; touch stamp/$c.Stop)\""
	done > jobs

	$pythonExec $jobsExec jobs

	echo Finished Pseudopipe on strand $t

	cd ..
done


# (7) Generating final results
outFilePrefix=$outDir/pgenes/`basename $outDir`_pgenes
$genPgeneResult $outDir $outFilePrefix.txt
awk '!a[$1$2$3$4$5]++' ${pre}_pgenes.txt  > ${pre}_pgenes_final.txt
$genFullAln     $outDir $outFilePrefix.align.gz

