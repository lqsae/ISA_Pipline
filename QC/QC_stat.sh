#!/bin/bash
if [ $# != 4 ]
then
	echo ""
	echo "  Version:         v2.01 (2017-12-18)"
	echo "  Authors:         Haibo Li"
	echo "  Description:     stat Fastq data basic infomation"
	echo ""
	echo "  Usage:           bash $0 <R1.fastq.gz list | R1.fastq list> <R2.fastq.gz list | R2.fastq list>  <Output_dir>  <SampleName> [thread]"
	echo ""
	echo "  <R1.fastq.gz list | R1.fastq list> :	R1 fastq[.gz] file list, seperated by comma."
	echo "  <R2.fastq.gz list | R2.fastq list> :	R2 fastq[.gz] file list, seperated by comma."
	echo "  [thread] :    specify cpu thread count ,default: 6."
	echo ""
	exit 1
fi

R1=$1
R2=$2
Output_dir=$3
work_dir=$Output_dir
Sample=$4

Step_Name="FastqStat"

#source software environment
Pipeline=$(cd `dirname $0`; pwd)

source $Pipeline/software.sh
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
nt=6

#-----------------------------------------------------------------------
# 0. Check fastq file
#-----------------------------------------------------------------------
# 0.1 check R1 file & count
#-----------------------------------------------------------------------
log_event 0.1 "check R1 file" >> $step_log
R1_list=${R1//,/ }
R1_fq_list=""
R1_count=0
for R1_file in ${R1//,/ }
do
	let R1_count=$R1_count+1
	if [ ! -f $R1_file ]
	then
		echo "Error: $R1_file not exists"
		exit 1;
	fi
	R1_fq_list=$R1_fq_list"$R1_file "
done
#-----------------------------------------------------------------------
# 0.2 check R2 file & count
#-----------------------------------------------------------------------
log_event 0.2 "check R2 file" >> $step_log
R2_list=${R2//,/ }
R2_fq_list=""
R2_count=0
for R2_file in ${R2//,/ }
do
	let R2_count=$R2_count+1
	if [ ! -f $R2_file ]
	then
		echo "Error: $R2_file not exists"
		exit 1;
	fi
	R2_fq_list=$R2_fq_list"$R2_file "
done
#-----------------------------------------------------------------------
# 0.3 check R1 count  equal R2 count
#-----------------------------------------------------------------------
log_event 0.3 "check R1 & R2 equal Number" >> $step_log
if [ $R1_count != $R2_count ]
then
	echo "Error: R1 fastq file count not equal R2 fastq file count"
	exit 1;
fi

#-----------------------------------------------------------------------
# 0.4 merge R1 or R2 file
#-----------------------------------------------------------------------
log_event 0.4 "Merge multiple R1 or R2 fastq" >> $step_log
R1_name=""
R2_name=""
R1_format=""
R2_format=""
R1_raw=""
R2_raw=""
if [ $R1_count == 1 ]
then
        R1_name=`basename $R1_list | sed 's/.fastq.*//g;s/.fq.*//g'`
	R1_format=`basename $R1_list | sed 's/.*.fastq//g;s/.*.fq//g'`
	ln -s $R1_list $Sample.raw.R1.fastq$R1_format
	R1_raw=$Sample.raw.R1.fastq$R1_format

        R2_name=`basename $R2_list | sed 's/.fastq.*//g;s/.fq.*//g'`
        R2_format=`basename $R2_list | sed 's/.*.fastq//g;s/.*.fq//g'`
	ln -s $R2_list $Sample.raw.R2.fastq$R2_format
	R2_raw=$Sample.raw.R2.fastq$R2_format
else
	tmp_R1_file=`echo $R1_list | sed 's/ .*//g'`
	R1_format=`basename $tmp_R1_file | sed 's/ .*//g;s/.*.fastq//g;s/.*.fq//g'`
	cat $R1_list > $Sample.raw.R1.fastq$R1_format
	R1_raw=$Sample.raw.R1.fastq$R1_format
	tmp_R2_file=`echo $R2_list | sed 's/ .*//g'`
	R2_format=`basename $tmp_R2_file | sed 's/.*.fastq//g;s/.*.fq//g'`
	cat $R2_list > $Sample.raw.R2.fastq$R2_format
	R2_raw=$Sample.raw.R2.fastq$R2_format
fi

#-----------------------------------------------------------------------
# Data statistics
#-----------------------------------------------------------------------
log_event qc_stat.2 "Data stat begin" >> $step_log

java -jar $qcstat/qcstat.jar 33 $Sample $R1_raw $R2_raw
perl $Pipeline/QC/stat_Qual.pl $Sample.R1.qual.stat $Sample.R2.qual.stat > $Sample.qual.stat
Rscript $Pipeline/QC/plot_Qual.R $Sample.qual.stat $Sample
Rscript $Pipeline/QC/plot_GC.R $Sample $Sample.R1 $Sample.R2

for i in $Sample*.pdf
do
	convert -density 100 -background white -flatten $i ${i%pdf}png
done

set +e
#-----------------------------------------------------------------------
# 1.4 Remove temp file
#-----------------------------------------------------------------------
log_event qc_stat.2 "remove temp file" >> $step_log
mkdir ${Sample}_temp
mv $R1_raw $R2_raw  ${Sample}_temp
mv ${Sample}*.{qmean,gc,qual,basic}.stat ${Sample}_temp

#-----------------------------------------------------------------------
log_event qc_stat.3 "$Sample $Step_Name finished" >> $step_log
