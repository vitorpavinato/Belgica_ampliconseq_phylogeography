#!/bin/bash

# Stop script on error
set -ue

# Configuration
bayescan="path_to_bayescan"
input_folder="path_to_input_folder"
input_file="input_file_name" 
remove_snps="list_of_snps_to_remove"
output_folder="path_to_output_folder"

./${bayescan} ${input_folder}/${input_file} -d ${input_folder}/${remove_snps} -od ${output_folder} -o newrun -all_trace -threads 30 -n 500000 -nbp 20 -pilot 50000 -burn 500000 -out_pilot -out_freq

# remove SNPs present in the mtDNA
#bayescan_2.1 belgica_1260_snps_bayescan.txt -d remove_bayescan_mtDNA.txt \ 
#                                            -od /home/pavinato/Dropbox/PosDoc_OSU_2015/Belgica_antarctica_genomics/WGS_analysis_pipeline_1/results/detect_selection/outlier_detection/bayescan/run \ 
#                                            -o newrun -all_trace -threads 30 -n 500000 -nbp 20 -pilot 50000 -burn 500000 -out_pilot -out_freq
