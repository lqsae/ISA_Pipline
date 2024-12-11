#!/bin/perl
if(@ARGV!=2){
	print  "\n Usage:  perl $0 <SampleName list> <Sample.dedup_metrics.txt list>\n";
	print "\n\tNote:\n";
	print "\t\tlist must be seperated by comma ','\n\n";
	exit 1;
}

my @sample = split(",",$ARGV[0]);
my @file = split(",",$ARGV[1]);
if(@sample != @file){ die "\nError: SampleName list count != Sample.dedup_metricx.txt file list\n\n";}

print "Sample\tTotalReads\tMapped Reads\tMapRatio(%)\tDupReads\tDupRatio(%)\n";
for(my $i=0;$i<@sample;$i++){
	my $name=$sample[$i];
	open FILE,"$file[$i]";
	my $header="";
	while(<FILE>){
		next if /^#/;
		chomp;
		if($_=~/^LIBRARY/){
			$header=$_;
			last;
		}
	}
	my %index;
	my @first_line=split("\t",$header);
	for(my $col=0;$col<@first_line;$col+=1){
		$index{$first_line[$col]}=$col;
	}
	my $value=<FILE>;
	chomp $value;
	my @second_line=split("\t",$value);
	my $mapped_reads = $second_line[$index{"READ_PAIRS_EXAMINED"}] * 2 + $second_line[$index{"UNPAIRED_READS_EXAMINED"}];
	my $unmap_reads = $second_line[$index{"UNMAPPED_READS"}];
	my $dup_reads = $second_line[$index{"READ_PAIR_DUPLICATES"}] * 2 + $second_line[$index{"UNPAIRED_READ_DUPLICATES"}];
	my $total_reads=$mapped_reads + $unmap_reads;
	printf "%s\t%d\t%d\t%.3f\t%d\t%.3f\n", $name,$total_reads,$mapped_reads,$mapped_reads*100.0/$total_reads,$dup_reads,$dup_reads*100.0/$mapped_reads;
	close FILE;

}


