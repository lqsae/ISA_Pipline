#!/bin/perl
use List::Util qw/sum/;
if(@ARGV!=4){
	print "\n  Usage:  perl $0 <Sample.Depth.freq>  <RefGenome size>  <Depth_dist list>  <Sample_name>\n\n";
	exit 1;
}

my $prefix=$ARGV[3];
my $depth_dist=$ARGV[2];
my $total_size=$ARGV[1];

open DP,"$ARGV[0]";

my @depth_freq;
my $cov_size=0;
my $sum_depth=0;
while(<DP>){
	next if /^Depth/;
	chomp;
	my @tmp_line=split("\t",$_);
	$depth_freq[$tmp_line[0]]=$tmp_line[1];
	$sum_depth += ($tmp_line[0] * $tmp_line[1]);
	if($tmp_line[0]>0){$cov_size+=$tmp_line[1];}
}
close DP;

my @depth_list=split(",",$depth_dist);

my $header_info="Sample\tAve_Depth(X)";
my $summary_info="$prefix\t".sprintf("%.2f",$sum_depth/$cov_size);

for(my $i=0;$i<@depth_list;$i++){
	$header_info.="\t>=".$depth_list[$i]."X";
	my $tmp_size=sum(@depth_freq[$depth_list[$i]..$#depth_freq]);
	$summary_info.="\t".sprintf("%.2f",$tmp_size*100.0/$total_size)."%";
}

print $header_info."\n";
print $summary_info."\n";
