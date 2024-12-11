#!/bin/bash

#----------模块信息------------------
#  Version:	v1.02 (2022-01-21)"
#  Authors:	Haibo Li
#------------------------------------
# 模块帮助文档
function usage(){
	echo ""
	echo "  Description:     对Pairend模式Fastq文件进行质控，去除接头，去除低质量，去除短序列"
	echo ""
	echo "  Usage:	bash $0 [options]  <R1.fastq.gz list|R1.fastq list>  <R2.fastq.gz list|R2.fastq list>  <Read1 Adapter_seq>  <Read2 Adapter_seq>  <输出目录>  <输出前缀>"
	echo ""
	echo "		<R1.fastq.gz list|R1.fastq list> :	R1 fastq[.gz] file list, seperated by comma."
	echo "		<R2.fastq.gz list|R2.fastq list>> :	R2 fastq[.gz] file list, seperated by comma."
	echo "		<Read1 Adapter_seq> :   		R1接头序列"
	echo "		<Read2 Adapter_seq> :			R2接头序列"
	echo "		<输出目录> :				指定输出目录"
	echo "		<输出前缀> :				指定输出前缀，建议样本编号，样本_trio，家系编号_trio"
	echo ""
	echo "  Options:  选项参数说明"
	echo ""
	echo "		-t :			可选参数，指定线程数, 默认值：6"
	echo ""
	echo "		--stepname :		自定义任务名, e.g. QC, Align, SNV....， 用于区分不同模块的日志文件和结束标识文件，默认值：脚本的前缀"
	echo "		--pipeline :		指定本模块所在路径，默认值：脚本真实路径目录. 暂未启用"
	echo "		--module :		指定本模块依赖的脚本目录，默认值：脚本真实路径目录/脚本前缀"
	echo "		--tmp :			指定该选项时，则保留中间临时文件_temp目录"
	echo ""
	echo "  Note:	备注说明信息"
	echo "	    <Adapter_seq> default:"
	echo "		Illumina:     AGATCGGAAGAGC AGATCGGAAGAGC"
	echo "		Nextera:      CTGTCTCTTATA CTGTCTCTTATA"
	echo "		BGI:          AAGTCGGAGGCCA AAGTCGGATCGTA"
	echo ""
	echo "  software default: cutadapt"
        echo ""
	exit 1
}

# 参数设置， 需自定义	
ARGS=`getopt -o t:h --long stepname:,pipeline:,module:,tmp -n "$0" -- "$@"`

if [ $? != 0 ] || [ $# == 0 ] ; then
	usage
fi
eval set -- "${ARGS}"
 
# 设置变量默认值, 需自定义
nt=6

Step_Name="`basename $(realpath $0) | sed 's/\..*//g'`"		# 初始化任务名
Pipeline=$(dirname $(realpath $0))				# 初始化流程/模块路径
Module=$Pipeline/$Step_Name					# 初始化流程/模块依赖脚本目录名，需存在于Pipeline目录中

Clean_tmp="Y"

# 选项参数传递，与参数核查， 需自定义
while true
do
    case "$1" in
	-h)		usage ;			exit 0	;;
	-t)		nt=$2 ;			shift 2	;; 
	--stepname)	Step_Name=$2;		shift 2	;;
	--pipeline)	Pipeline=$(realpath $2);shift 2	;;
	--module)	Module=$(realpath $2);	shift 2 ;;
	--tmp)		Clean_tmp="N";		shift 1 ;;
        --)		shift ;			break	;;
	?)       	echo "不明参数 $1" ;    exit 1	;;
        *)		echo "内部错误!" ;	exit 1	;;
    esac
done

# 选项参数确认
if [ ! -d $Module ] ; then echo "Error: 未找到依赖的脚本目录 $Module, 请指定--module参数"; exit 1; fi

# 位置参数传递， 需自定义
if [ $# != 6 ] ; then echo "Error: 位置参数不符！" ; usage ; fi

R1=$1
R2=$2
adapter1=$3
adapter2=$4

Output_dir=$5
work_dir=$Output_dir
Sample=$6

# 模块运行开始
set -e
Run_begin=`date +%s`

# 导入日志函数
source $Module/log.sh
# 导入软件环境
source $Module/software.sh

# 创建输出目录，并定位该目录
mkdir -p $work_dir
cd $work_dir

# 创建日志
step_log=$work_dir/$Sample.$Step_Name.steps.log
log_Command $0 $(echo $ARGS|sed "s/ -- / /g;s/'//g") > $step_log

detail_log=$work_dir/$Sample.$Step_Name.details.log
exec >${detail_log} 2>&1

set -x

# 位置参数确认, 自定义
check_file "R1_FASTQ" $R1
check_file "R2_FASTQ" $R2

# 导入数据库环境
#source $Module/hg38.env.sh

#-----------------------------------------------------------------------
# 方法/函数主体, 自定义
#-----------------------------------------------------------------------
# 第1步:	确认R1 & R2
#-----------------------------------------------------------------------
log_event auto.1 "Check fastq begin" >> $step_log

R1_list=${R1//,/ }
R1_fq_list=""
R1_count=0
for R1_file in ${R1//,/ }
do
	let R1_count=$R1_count+1
	R1_fq_list=$R1_fq_list"$R1_file "
done

#----------------------------------
R2_list=${R2//,/ }
R2_fq_list=""
R2_count=0
for R2_file in ${R2//,/ }
do
	let R2_count=$R2_count+1
	R2_fq_list=$R2_fq_list"$R2_file "
done

#----------------------------------
if [ $R1_count != $R2_count ]
then
	echo "Error: R1 fastq file count not equal R2 fastq file count"
	exit 1;
fi

#-----------------------------------------------------------------------
# 第2步:	文件合并
#-----------------------------------------------------------------------
log_event auto.2 "Merge fastq begin" >> $step_log

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
	ln -f -s $R1_list $Sample.raw.R1.fastq$R1_format
	R1_raw=$Sample.raw.R1.fastq$R1_format
	R2_name=`basename $R2_list | sed 's/.fastq.*//g;s/.fq.*//g'`
	R2_format=`basename $R2_list | sed 's/.*.fastq//g;s/.*.fq//g'`
	ln -f -s $R2_list $Sample.raw.R2.fastq$R2_format
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
# 第3步:     Fastq文件质控
#-----------------------------------------------------------------------
log_event auto.3 "Trim fastq begin" >> $step_log

cutadapt -a $adapter1 -A $adapter2 -e 0.1 -O 1 -m 36 --max-n 5 -q 30,30 -j $nt -o $Sample.trim.R1.fastq.gz -p $Sample.trim.R2.fastq.gz $R1_raw $R2_raw > $Sample.trim.log

#-----------------------------------------------------------------------
# 第4步:     Fastq统计
#-----------------------------------------------------------------------
log_event auto.4 "Stat fastq begin" >> $step_log

java -jar $qcstat/qcstat.jar 33 $Sample.Raw $R1_raw $R2_raw &
Raw_QC_pid=$!
java -jar $qcstat/qcstat.jar 33 $Sample.Clean $Sample.trim.R1.fastq.gz $Sample.trim.R2.fastq.gz &
Clean_QC_pid=$!
wait $Raw_QC_pid $Clean_QC_pid

perl $Module/stat_Qual.pl $Sample.raw.R1.qual.stat $Sample.raw.R2.qual.stat > $Sample.raw.qual.stat
perl $Module/stat_Qual.pl $Sample.trim.R1.qual.stat $Sample.trim.R2.qual.stat > $Sample.trim.qual.stat
# Rscript $Module/plot_Qual.R $Sample.raw.qual.stat $Sample.Raw
# Rscript $Module/plot_Qual.R $Sample.trim.qual.stat $Sample.Clean
# Rscript $Module/plot_GC.R $Sample.Raw $Sample.raw.R1 $Sample.raw.R2
# Rscript $Module/plot_GC.R $Sample.Clean $Sample.trim.R1 $Sample.trim.R2

for i in $Sample*.pdf
do
	convert -density 100 -background white -flatten $i ${i%pdf}png
done

cat $Sample.Raw.stat $Sample.Clean.stat | grep -v "Reads" | awk -F'\t' -v sample=$Sample 'BEGIN{print "Sample\tRawReads\tCleanReads\tReadsRatio(%)\tRawBases\tCleanBases\tBasesRatio(%)"}; \
        {if(NR==1){raw_R=$2;raw_B=$3} \
        if(NR==2){clean_R=$2;clean_B=$3}} \
        END{R_ratio=clean_R*100.0/raw_R;B_ratio=clean_B*100.0/raw_B; \
        printf "%s\t%d\t%d\t%.2f\t%d\t%d\t%.2f\n",sample,raw_R,clean_R,R_ratio,raw_B,clean_B,B_ratio}' > $Sample.Ratio.stat

set +e
#-----------------------------------------------------------------------
# Remove temp file
#-----------------------------------------------------------------------
log_event auto.clean "Clean temp file" >> $step_log
mkdir ${Sample}_temp
mv $R1_raw $R2_raw  ${Sample}_temp
mv ${Sample}*trim.log ${Sample}_temp
mv ${Sample}*.{qmean,gc,qual,basic}.stat ${Sample}_temp

if [ $Clean_tmp != "N" ]
then
	rm -rf ${Sample}_temp
fi

#-----------------------------------------------------------------------
log_event auto.final "$Sample $Step_Name finished" >> $step_log

set +x 
Run_end=`date +%s`
check_times "$Sample $Step_Name Run Times:" $Run_begin $Run_end > $work_dir/$Sample.$Step_Name.done
