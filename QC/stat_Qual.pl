#!/bin/perl
use Statistics::Descriptive;
if(@ARGV!=2){
	print "\n  Usage:  perl $0 <R1.qual list> <R2.qual list>\n\n";
	exit 1;
}

my @R1_array=split(",",$ARGV[0]);
my @R2_array=split(",",$ARGV[1]);

my @R1_qual_array;
my @R2_qual_array;

foreach my $file (@R1_array){
	open IN,"$file";
	while(<IN>){
		chomp;
		my @tmp_line=split("\t",$_);
		my $pos=$tmp_line[0];
		for(my $i=1;$i<@tmp_line;$i++){
			my $qual=$i+1;
			my $count=$tmp_line[$i];
			$R1_qual_array[$pos][$qual]+=$count;
		}
	}
	close IN;
}

foreach my $file (@R2_array){
	open IN,"$file";
	while(<IN>){
		chomp;
		my @tmp_line=split("\t",$_);
		my $pos=$tmp_line[0];
		for(my $i=1;$i<@tmp_line;$i++){
			my $qual=$i+1;
			my $count=$tmp_line[$i];
			$R2_qual_array[$pos][$qual]+=$count;
		}
	}
	close IN;
}

print "Pos\tR1_ave_qual\tR1_ave_error\tR2_ave_qual\tR2_ave_error\n";
my $max_pos=(@R1_qual_array>@R2_aual_array)?@R1_qual_array:@R2_qual_array;

for(my $pos=1;$pos<$max_pos;$pos++){
	my $R1_sum_count=0;
	my $R1_sum_qual=0;
	my $R2_sum_count=0;
	my $R2_sum_qual=0;
	for(my $qual=2;$qual<=41;$qual++){
		my $R1_count=$R1_qual_array[$pos][$qual];
		my $R2_count=$R2_qual_array[$pos][$qual];

		$R1_sum_count+=$R1_count;
		$R2_sum_count+=$R2_count;
		
		$R1_sum_qual+=$R1_count*$qual;
		$R2_sum_qual+=$R2_count*$qual;

#		print $pos."\t".$qual."\t".$R1_sum_count."\t".$R1_sum_qual."\n";
	}
#	print $R1_sum_count."\t".$R1_sum_qual."\n";
	my $R1_ave_qual=$R1_sum_qual*1.0/$R1_sum_count;
	my $R2_ave_qual=$R2_sum_qual*1.0/$R2_sum_count;
	my $R1_ave_error=10**($R1_ave_qual/10*-1);
	my $R2_ave_error=10**($R2_ave_qual/10*-1);
	printf "%d\t%e\t%e\t%e\t%e\n",$pos,$R1_ave_qual,$R1_ave_error,$R2_ave_qual,$R2_ave_error;
}

