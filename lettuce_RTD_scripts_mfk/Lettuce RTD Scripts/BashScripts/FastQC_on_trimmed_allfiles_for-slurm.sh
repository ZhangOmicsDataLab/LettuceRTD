#! /bin/bash

for fq in ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/results/03_fastq_trimmed/*_trimmed.fq.gz
do

sbatch -p nodes -t 00:15:00 -c 1 --job-name fastqc_trimmed --mem 15G --output ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/logs/output.%j.out --wrap="sh ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/scripts/FastQC_on_trimmed_input_file.sh $fq"

sleep 1	# wait 1 second between each job submission
  
done
