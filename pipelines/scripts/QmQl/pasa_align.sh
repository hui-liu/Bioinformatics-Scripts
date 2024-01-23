cp ~/bin/PASApipeline-v2.5.2/pasa_conf/pasa.alignAssembly.Template.txt alignAssembly.config
cp ~/bin/PASApipeline-v2.5.2/pasa_conf/pasa.annotationCompare.Template.txt annotationCompare.config
mkdir db

dir=$1
sed -i "s#<__DATABASE__>#$dir/db/pasa.sqlite#" alignAssembly.config
sed -i "s#<__DATABASE__>#$dir/db/pasa.sqlite#" annotationCompare.config

ln -s $2 $3.fasta

echo "~/bin/PASApipeline-v2.5.2/Launch_PASA_pipeline.pl \
-c alignAssembly.config -R -g $3.fasta \
-t $3_transcriptome.fasta.clean \
-T -u $3_transcriptome.fasta \
--ALIGNERS gmap,blat \
--CPU 48 --stringent_alignment_overlap 30.0 \
--MAX_INTRON_LENGTH 60000 --TRANSDECODER -C" > run_alignAssembly.sh
echo "cp -r db db.raw" >> run_alignAssembly.sh
