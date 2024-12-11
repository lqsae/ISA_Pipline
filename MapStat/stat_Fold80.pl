#!/bin/perl
if(@ARGV!=3){
	print "\n  Usage:  perl $0 <Depth.freq file> <RefGenome size>  <Sample_name>\n\n";
	exit 1;
}
my $prefix=$ARGV[2];
my $total_size=$ARGV[1];

open DP,"$ARGV[0]";
my @depth_freq;
my $cov_size=0;
my $sum_depth=0;
while(<DP>){
	chomp;
	my @tmp_line=split("\t",$_);
	$depth_freq[$tmp_line[0]] = $tmp_line[1];
	$sum_depth += $tmp_line[0] * $tmp_line[1];
	if($tmp_line[0]>0){$cov_size += $tmp_line[1];}
}
close DP;

my $cumul_size=0;
my ($depth10_cov,$depth20_cov,$depth30_cov,$depth50_cov)=(0,0,0,0);
my $dp_percent80="";
for(my $i=1;$i<@depth_freq;$i+=1){
	if($depth_freq[$i] eq ""){$depth_freq[$i]=0;}
	$cumul_size+=$depth_freq[$i];
	my $cumul_ratio=sprintf("%.3f",$cumul_size*100.0/$cov_size);
	print DIST $i."\t".$depth_freq[$i]."\t".sprintf("%.3f",$depth_freq[$i]*100.0/$cov_size)."\t".$cumul_ratio."\n";
	if($dp_percent80 eq "" && $cumul_ratio >= 20){$dp_percent80 = $i; }
	if($i==9){
		$depth10_cov= sprintf("%.2f",(($cov_size-$cumul_size)*100.0/$total_size));
	}
	if($i==19){
		$depth20_cov= sprintf("%.2f",(($cov_size-$cumul_size)*100.0/$total_size));
	}
	if($i==29){
		$depth30_cov= sprintf("%.2f",(($cov_size-$cumul_size)*100.0/$total_size));
	}
	if($i==49){
		$depth50_cov= sprintf("%.2f",(($cov_size-$cumul_size)*100.0/$total_size));
	}

}

print "Sample\tAve_Depth(X)\tFold80\t>=1X\t>=10X\t>=20X\t>=30X\t>=50X\n";
print "$prefix\t".sprintf("%.2f",$sum_depth/$cov_size)."\t".sprintf("%.3f",($sum_depth/$cov_size)/$dp_percent80)."\t".sprintf("%.2f",$cov_size*100.0/$total_size)."%\t".$depth10_cov."%\t".$depth20_cov."%\t".$depth30_cov."%\t".$depth50_cov."%\n";
close COV;
close DIST;

