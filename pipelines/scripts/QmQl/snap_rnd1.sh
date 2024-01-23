maker2zff -x 0.25 -l 50 -d ../../maker/rnd1/master_datastore_index.log

# gather some stats and validate
fathom genome.ann genome.dna -gene-stats > gene-stats.log
fathom genome.ann genome.dna -validate > validate.log

# collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom genome.ann genome.dna -categorize 1000 > categorize.log
fathom uni.ann uni.dna -export 1000 -plus > export.log

# create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log
cd ..

# assembly the HMM
hmm-assembler.pl genome params > genome.hmm

