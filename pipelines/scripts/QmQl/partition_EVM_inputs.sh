dir=$1
genome=$2

awk '$2=="maker"' $dir/maker/rnd3/rnd3_noSeq.gff > maker.gff

python2 ~/bin/scripts/maker2gff3.py \
$dir/maker/rnd3/rnd3_noSeq.gff \
augustus augustus.gff

python2 ~/bin/scripts/maker2gff3.py \
$dir/maker/rnd3/rnd3_noSeq.gff \
snap snap.gff

python2 ~/bin/scripts/maker2gff3.py \
$dir/maker/rnd3/rnd3_noSeq.gff \
genemark genemark.gff

cat maker.gff \
augustus.gff \
snap.gff \
genemark.gff > gene_predictions.gff3

awk '$3=="match_part"' $dir/maker/rnd3/rnd3.protein2genome.gff | \
sed 's/ID=.*Parent/ID/' | sed 's/match_part/protein_match/' > protein_alignments.gff3

awk '$3=="match_part"' $dir/maker/rnd3/rnd3.est2genome.gff | \
sed 's/ID=.*Parent/ID/' | \
sed 's/match_part/expressed_sequence_match/' > transcript_alignments.gff3

# check
~/bin/PASApipeline-v2.5.2/misc_utilities/pasa_gff3_validator.pl gene_predictions.gff3
~/bin/PASApipeline-v2.5.2/misc_utilities/pasa_gff3_validator.pl protein_alignments.gff3
~/bin/PASApipeline-v2.5.2/misc_utilities/pasa_gff3_validator.pl transcript_alignments.gff3


# weights
echo -e "TRANSCRIPT\test2genome\t10
PROTEIN\tprotein2genome\t5
ABINITIO_PREDICTION\tmaker\t4
ABINITIO_PREDICTION\taugustus\t2
ABINITIO_PREDICTION\tsnap\t2
ABINITIO_PREDICTION\tgenemark\t1" > weights.txt


echo "partition_EVM_inputs.pl --segmentSize 1000000 --overlapSize 200000 \
--genome $genome \
--gene_predictions gene_predictions.gff3 \
--protein_alignments protein_alignments.gff3 \
--transcript_alignment transcript_alignments.gff3 \
--partition_listing partitions_list.out

write_EVM_commands.pl \
--genome $genome \
--gene_predictions gene_predictions.gff3 \
--protein_alignments protein_alignments.gff3 \
--transcript_alignment transcript_alignments.gff3 \
--weights $dir/evm/weights.txt \
--partitions partitions_list.out \
--output_file_name evm.out > EVM_commands.list" > run_partition_EVM_inputs.sh
