#! /usr/bin/bash

#if [$# -eq 0];then
	echo "Usage:bash run_dada.sh -i sample_data_path -o ps_file_name -r R_data"
	echo "eg: bash run_dada.sh -i /home/jianglikun/dada2/sum_data/ -o ps -r R.data"
#	exit 1
#fi

while getopts ":i:o:r:h" arg
do
	case $arg in
		i)
			sample_data=$OPTARG;;
		o)
			output_ps=$OPTARG;;
		r)
			R_data=$OPTARG;;
		h)
			echo "Usage:run_dada.sh -i sample_dada -o ps_file_name -r R_data "
			echo "eg: bash run_dada.sh -i /home/jianglikun/dada2/sum_data/ -o ps -r test_R_data"
			exit 1;;
	esac
done

/home/jianglikun/Toolkits/R-3.4.0/bin/Rscript run_dada.R $sample_data $outputfile $R_data
