#!/bin/bash
if [ $# != 4 ] && [ $# != 5 ]
then
        echo ""
        echo "  Version:         v2.01 (2017-11-28)"
        echo "  Authors:         Haibo Li"
        echo "  Description:     statistic BAM file basic information"
        echo ""
        echo "  Usage:           bash $0  <Sample.bam file>  <hg19|hg38>  <Output_dir>  <SampleName> [Pannel_BED file]"
        echo ""
        echo ""
        exit 1
fi

BAM_file=$1
Ref_version=$2
Output_dir=$3
work_dir=$Output_dir
Sample=$4
Bed_file=$5

Step_Name="AlignStat"

#source software environment
Pipeline=$(cd `dirname $0`; pwd)

source $Pipeline/Source/software.sh
source $Pipeline/Source/log.sh

# creat work directory
mkdir -p $work_dir
cd $work_dir

# creat log file
step_log=$work_dir/$Sample.$Step_Name.steps.log
log_Command $0 $@ > $step_log

detail_log=$work_dir/$Sample.$Step_Name.details.log
exec >${detail_log} 2>&1
set -e -x

# cpu number used
nt=4

# detect Bed file exists
PanelSize=""
if [ $Bed_file ]
then
	check_file "BED" $Bed_file
	PanelSize=`awk -F'\t' '{sum+=($3-$2)}END{print sum}' $Bed_file`
fi

# reference databases

if [ $Ref_version != "hg19" ] && [ $Ref_version != "hg38" ]
then
	echo "Error: $Ref_version version not correct, please specified hg19 or hg38"
	exit 1;
fi

source $Pipeline/Source/$Ref_version.env.sh
GenomeSize=`awk -F'\t' '{sum+=$2;}END{print sum}' $Ref_fai`

#-----------------------------------------------------------------------
# 2. Mapping and sorting
#-----------------------------------------------------------------------

#Off_pid=""
#Cov_pid=""
#Depth_pid=""

if [ $Bed_file ]
then
	
	nohup bash $Pipeline/MapStat/OffStat.sh $BAM_file $Bed_file $Ref_version $work_dir $Sample &
        Off_pid=$!

	nohup python $Pipeline/MapStat/coverage.py $BAM_file $Bed_file $RefGenome $Sample Y &
	Cov_pid=$!

	nohup samtools depth -d 100000 $BAM_file -b $Bed_file > ${Sample}.Panel.depth &
	Depth_pid=$!


	wait $Depth_pid
        perl $Pipeline/MapStat/stat_Cov.pl ${Sample}.Panel.depth $PanelSize $Sample.Panel
	Rscript $Pipeline/MapStat/plot_Depth-Dis.R $Sample.Panel.Depth.freq $Sample
else
	samtools depth -d 100000 $BAM_file  > ${Sample}.Genome.depth
	perl $Pipeline/MapStat/stat_Cov.pl ${Sample}.Genome.depth $GenomeSize $Sample.Genome
	Rscript $Pipeline/MapStat/plot_Depth-Dis.R $Sample.Genome.Depth.freq $Sample
fi

#-----------------------------------------------------------------------
# 2.6 infer gender
#-----------------------------------------------------------------------
log_event 2.6 "Infer gender begin" >> $step_log
samtools depth -r $SRY_Region $BAM_file > $Sample.$Ref_version.SRY.depth
filesizeb=`awk 'BEGIN{sum=0;}{sum+=$3};END {print sum }' $Sample.$Ref_version.SRY.depth `

#threldhold=5000
echo $filesizeb
Sex=""
if [ "$filesizeb" -lt 5000 ]
then
	Sex="female"
else
	Sex="male"
fi
echo -e "$Sample\t$Sex" > $Sample.sexinfer.stat

wait $Off_pid $Cov_pid $Depth_pid

set +e

#-----------------------------------------------------------------------
# 2.8 remove temp file
#-----------------------------------------------------------------------
log_event 2.7 "clean temp begin" >> $step_log

mkdir ${Sample}_temp
#rm ${Sample}*.depth

log_event 2.8 "$Sample alignment finished" >> $step_log
