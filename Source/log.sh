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
		echo $tmp_file
		if [ ! -f $tmp_file ]
		then
			echo "Error: $type file '$tmp_file' not exists"
			exit 1;
		fi
	done
}

function check_faminfo(){
	local fam_info=$1

	if [ ! -f $fam_info ]
	then
	        echo "Error: Family info file '$fam_info' not exists"
	        exit 1
	fi

        local col_count=`awk -F'\t' '{print NF}' $fam_info|sort|uniq`
        if [ $col_count != 6 ] && [ $col_count != 7 ]
        then
       	        echo "Error: Family info file '$fam_info' format not correct, must have 6|7 columns each line"
		echo -e "\t\t#test_id\tsample_id\tgroup\tgender\tR1_list\tR2_list"
               	exit 1
        fi

	if [ $col_count == 6 ]
	then
		while read test_id sample_id group gender R1_list R2_list
		do
			if [ "$gender" != "XX" ] && [ "$gender" != "XY" ] && [ "$gender" != "Male" ] && [ "$gender" != "Female" ] && [ "$gender" != "gender" ] && [ "$gender" != "-" ]
			then
				echo "Error: "$fam_info" file $sample_id info '$gender'  error, must be [XX|XY|Male|Female|-] "
				exit 1
			fi
		done < $fam_info 
	fi
}
