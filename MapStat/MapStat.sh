#!/bin/bash
if [ $# != 4 ] && [ $# != 5 ]
then
        echo ""
        echo "  Version:         v2.01 (2017-11-28)"
        echo "  Authors:         Haibo Li"
        echo "  Description:     Statistic Coverage and Depth "
        echo ""
        echo "  Usage:           bash $0  <Sample.realigned.bam> <RefGenome fai file> <Output_dir>  <SampleName> [Pannel_BED file]"
        echo ""
        echo ""
        exit 1
fi


Bam_file=$1
Ref_fai=$2
Output_dir=$3
work_dir=$Output_dir
Sample=$4
Bed_file=$5

Step_Name="StatCov"
# Directory for pipeline existing
#Pipeline=$(cd `dirname $0`; pwd)
Pipeline=/work/Pipeline/DNA_WES/V2

#source software environment
source $Pipeline/software.sh

# detect Bed file exists
if [ $Bed_file ]
then
	if [ ! -f $Bed_file ]
	then
		echo "$Bed_file not exists!"
		exit 1;
	fi
fi

GenomeSize=`awk -F'\t' '{sum+=$2;}END{print sum}' $Ref_fai`

# cpu number used
nt=10

mkdir -p $work_dir
cd $work_dir

#-----------------------------------------------------------------------------------------------
#  Create log file
#-----------------------------------------------------------------------------------------------
# absolute file name of step log of project
step_log=$work_dir/$Sample.$Step_Name.steps.log
echo "************************************" > $step_log
echo "Command:    $0 $@" >> $step_log
echo "************************************" >> $step_log
# absolute file name of detail log of project
detail_log=$work_dir/$Sample.$Step_Name.details.log

# function for creating a log entry
# @param1: step number
# @param2: event description
log_event(){
  local step=$1
  local desc="$2"
  local time_mod=$(date +%s)
  echo "`date`" >>  ${step_log}
  echo -e "\tTime: ${time_mod}\tStep: ${step}\tDesc: \"${desc}\"" >> ${step_log}
}

exec >${detail_log} 2>&1
set -x

#-----------------------------------------------------------------------
# 1.1 Call depth and statistic
#-----------------------------------------------------------------------
log_event 1.1 "Local realignment begin"
if [ $Bed_file ]
then
	PanelSize=`awk -F'\t' '{sum+=($3-$2)}END{print sum}' $Bed_file`
	samtools depth -d 100000 $Bam_file -b $Bed_file | perl $Pipeline/MapStat/stat_Cov.pl - $PanelSize $Sample.Panel
else
	samtools depth -d 100000 $Bam_file | perl $Pipeline/MapStat/stat_Cov.pl - $GenomeSize $Sample.Genome
fi

log_event 1.final "$Sample StatCov finished"











