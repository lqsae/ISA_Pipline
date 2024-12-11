#!/bin/bash

#-----------------------------------------------------------------------------------------------
#  Create Log Function
#-----------------------------------------------------------------------------------------------

# Record Command
log_Command(){
	arguments=$*
	echo "************************************"
	echo "Command:    $arguments"
	echo "************************************"
}

#  Create log file

# function for creating a log entry
# @param1: step number
# @param2: event description
log_event(){
	local step=$1
	local desc="$2"
	local time_mod=$(date +%s)
	echo "`date`"
	echo -e "\tTime: ${time_mod}\tStep: ${step}\tDesc: \"${desc}\""
}

function check_file(){
	local type=$1
	local file_list=$2
	for tmp_file in ${file_list//,/ }
	do
#		echo $tmp_file
		if [ ! -f $tmp_file ]
		then
			echo "Error: $type file '$tmp_file' not exists"
			exit 1;
		fi
	done
}

function check_dir(){
	local type=$1
	local dir_list=$2
	for tmp_dir in ${dir_list//,/ }
	do
#		echo $tmp_dir
		if [ ! -d $tmp_dir ]
		then
			echo "Error: $type directory '$tmp_dir' not exists"
			exit 1;
		fi
	done
}

function check_times(){
	local step=$1
	local time1="$2"
	local time2="$3"

	local use_times_min=`echo -e $time1"\t"$time2 | awk -F'\t' '{printf "%.2f\n", ($2-$1)/60}'`
	local use_times_hour=`echo -e $time1"\t"$time2 | awk -F'\t' '{printf "%.2f\n", ($2-$1)/3600}'`

	echo -e "${step}\t${use_times_hour} h\t${use_times_min} min"
}

function check_ref(){
	local ref_version=$1 
	if [ $ref_version != "hg19" ] && [ $ref_version != "hg38" ] && [ $ref_version != "RCRS" ]
	then
	        echo "Error: $ref_version version not correct, please specified hg19 or hg38 or RCRS"
        	exit 1;
	fi
}

function check_faminfo(){
	local fam_info=$1

        local col_count=`awk -F'\t' '{print NF}' $fam_info|sort|uniq`
        if [ $col_count -ne 7 ]
        then
       	        echo "Error: Family info file '$fam_info' format not correct, must have 7 columns each line"
		echo -e "\t\t#Family_ID\tRelation\tSample\tSick\tR1_list\tR2_list\tGender"
               	exit 1
        fi

	while read Fam_id Relation Sample Sick R1_list R2_list Gender
	do
		if [ "$Gender" != "Male" ] && [ "$Gender" != "Female" ] && [ "$Gender" != "Gender" ] && [ "$Gender" != "-" ]
		then
			echo "Error: "$fam_info" file $Sample info '$Gender'  error, must be [Male|Female|-] "
			exit 1
		fi
	done < $fam_info 
}

