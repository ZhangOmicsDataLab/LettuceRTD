#! /bin/bash

for fq in ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/raw_data/*_1.fastq.gz
do

sbatch -p nodes -t 00:45:00 -c 1 --job-name fastp --mem 15G --output ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/logs/output.%j.out --wrap="sh ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/scripts/fastp_on_input_file.sh $fq"

sleep 1	# wait 1 second between each job submission
  
done
