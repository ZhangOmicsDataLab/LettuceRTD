#!/bin/bash

for fq in ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/data/G2_1.fastq.gz
do


sbatch --job-name rnaseq-workflow --mail-type NONE --mail-user mfk513@york.ac.uk --ntasks 1 --cpus-per-task 12 --mem 6G --partition nodes --time 0-02:30 --output ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/logs/%j.log --account biol-gensac-2019 --wrap="sh ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/scripts/rnaseq_analysis_on_input_file_salmon_only.sh $fq ~/scratch/LsRTD_v1.0_valdn/references/LettuceRTDv1_GCA_000143535.4_ASM14353v4_merged/LettuceRTDv1_GCA_000143535.4_ASM14353v4_merged_salmon_index"


sleep 1 # wait 1 second between each job submission

done
