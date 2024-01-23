# Combining the Partitions
recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

# Convert to GFF3 Format
genome=$1
convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output_file_name evm.out --genome $genome

# merge the out files and gff files
cat */evm.out.gff3 > evm.all.gff3
gff3_file_to_proteins.pl evm.all.gff3 $genome prot > evm.PEP.fasta
gff3_file_to_proteins.pl evm.all.gff3 $genome CDS  > evm.CDS.fasta

