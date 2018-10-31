#!/bin/bash
# USAGE:
# run this script like as follow sample:
# > command_to_run.sh <mode> <path_to_source_file> <path_to_experiment_dir> <path_to_output_dir> <sub_experiment_separated_by_comma_wo_spacers> <reversed_barcode>
# 
# EXAMPLE:
# > ./command_to_run.sh /home/anton/backup/input/trip/RUN_2018-05-10/sample_S1_L001_R1_001.fastq.gz /home/anton/backup/input/trip/RUN_2018-05-10/run/Lib_25-32 /home/anton/backup/input/trip/RUN_2018-05-10/results control_e,control_m,control_n,normalization,expression,mapping rev
# 
MODE=$1
SOURCE=`echo $2 | sed -r "s/,/ /g"`
PYTHON='/usr/bin/python2.7'
TRIP='TripMain_0_2.py'
ANALYSIS_MAIN='/home/anton/data/TRIP/pyMPFA/AnalysisMain.py'
EXPERIMENT_DIR=$3
OUTPUT=$4
LABEL=`basename $EXPERIMENT_DIR`
SUB_EXPERIMENT=`echo $5 | sed -r "s/,/ /g"`
REVERSED_BC=$6
if [[ ${REVERSED_BC} -eq "rev" ]]; then
	REVERSED_BC="-rb"
else
	REVERSED_BC=""
fi
# Run execution of TRIP
for subs in $SUB_EXPERIMENT; do
	cd ${EXPERIMENT_DIR}/${subs}
	if [[ $subs = "expression" || $subs = "normalization" || $subs = "control_e" || $subs = "control_n" ]]; then
		options="-B"
	else
		options=""
	fi
	source_count=`echo "$SOURCE" | wc -w`
	if [ $source_count -eq 1 ]; then
		main_input="-i $SOURCE"
	elif [ $source_count -eq 2 ]; then
		main_input=`echo "$SOURCE" | sed -r "s/(.*) (.*)/-f \1 -r \2/g"`
	fi
	echo "Start pipe ${LABEL}_${subs}..."
	echo "$PYTHON $TRIP -m $MODE $main_input -o $OUTPUT -l ${LABEL}_${subs} ${REVERSED_BC} $options"
	$PYTHON $TRIP -m $MODE $main_input -o $OUTPUT -l ${LABEL}_${subs} ${REVERSED_BC} $options
	echo "...end pipe"
done
echo "End of all pipes"
echo "Start intermediate results analysis..."
# $PYTHON $ANALYSIS_MAIN
