sample=$1
echo "~/bin/PASApipeline-v2.5.2/Launch_PASA_pipeline.pl \
-c ../annotationCompare.config -A -g ../${sample}.fasta \
-t ../${sample}_transcriptome.fasta.clean \
--CPU 48 -L --annots ../../evm/evm.all.gff3" > run_annotationCompare_rnd1.sh
echo "cp -r ../db ../db.rnd1" >> run_annotationCompare_rnd1.sh

