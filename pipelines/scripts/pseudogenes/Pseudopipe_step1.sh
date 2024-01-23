# Prepare
pre=$1
genome=$2
gff=$3
pep=$4
dir=$5
# Input files configuration
cd $dir
mkdir -p {ppipe_output/$pre,ppipe_input/$pre}
mkdir -p ppipe_input/$pre/{dna,pep,mysql}
cp $2 ppipe_input/$pre/dna/
cp $4 ppipe_input/$pre/pep/

cd $dir/ppipe_input/$pre/dna/
# a list file for all unmasked dna divided into different chromosomes in FASTA format
python /mnt/Data_disk/liuhui/bin/scripts/pseudogenes/split_seqs_by_chr.py $genome $pre
cd $dir/ppipe_input/$pre/mysql/
python /mnt/Data_disk/liuhui/bin/scripts/pseudogenes/extract_exon_by_chr.py $gff


# Pseudopipe configuration
PSEUDOPIPE_HOME=/mnt/Data_disk/liuhui/bin/pgenes/pseudopipe
export pseudopipe=$PSEUDOPIPE_HOME/core/runScripts.py
export genPgeneResult=$PSEUDOPIPE_HOME/ext/genPgeneResult.sh
export genFullAln=$PSEUDOPIPE_HOME/ext/genFullAln.sh
export fastaSplitter=$PSEUDOPIPE_HOME/ext/splitFasta.py
export sqDummy=$PSEUDOPIPE_HOME/ext/sqDummy.py
export blastHandler=$PSEUDOPIPE_HOME/core/processBlastOutput.py
export extractExLoc=$PSEUDOPIPE_HOME/core/extractKPExonLocations.py

# Alignment tools configuration
export formatDB=/mnt/Data_disk/liuhui/bin/blast-2.2.25/bin/formatdb
export blastExec=/mnt/Data_disk/liuhui/bin/blast-2.2.25/bin/blastall
export fastaExec=/mnt/Data_disk/liuhui/bin/fasta34/tfasty34

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
$pythonExec $fastaSplitter $inputPEP $((($pepNum+359)/360)) 'split%04d'
for f in `ls split*` 
do
	echo "( cd $(pwd); touch stamps/${f}.Start ; ( $blastExec -p tblastn -m 8 -z 4.7e8 -e .1 -d $inputDNA -i $f -o output/${f}.Out ; touch stamps/${f}.Stamp ) >status/${f}.Status 2>&1 )" 
done > jobs
