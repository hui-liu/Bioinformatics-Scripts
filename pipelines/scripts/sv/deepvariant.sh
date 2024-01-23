export PATH='/media/nfs2/wbs/wbs2/bin/samtools1.10/bin':$PATH
singularity run \
~/bin/deepvariant_singularity/deepvariant_1.5.0.sif \
/opt/deepvariant/bin/run_deepvariant \
--model_type=PACBIO \
--ref=$1 \
--reads=$2 \
--sample_name $3 \
--output_vcf=$3.vcf.gz \
--output_gvcf=$3.g.vcf.gz \
--intermediate_results_dir intermediate_results_dir \
--num_shards=96

rm intermediate_results_dir/make_examples.tfrecord-000*gz

