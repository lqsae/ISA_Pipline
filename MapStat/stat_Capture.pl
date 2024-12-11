#!/bin/perl
if(@ARGV!=3){
	print  "\n Usage:  perl $0 <Target.bed> <Out-Target.bed> <Read_Align.bed>\n";
	print "\n  Note:\n";
	print "        <Read_Align.bed>:\tcreated by bamToBed command.\n";
	print "        \t\t\tcommand: 'bamToBed -i Sample.bam > Sample.Reads.bed'\n\n";
	exit 1;
}

open BED,"$ARGV[0]";
open OUT_BED,"$ARGV[1]";
open READ_pos,"$ARGV[2]";

#system("date");
my %in_hash;
my %out_hash;
while(<BED>){
    my ($chr,$start,$end,$info) = split("\t");
    my $bin1 = int($start/10000);
	my $bin2 = int($end/10000);
	foreach my $k($bin1..$bin2){
		$in_hash{$chr}{$k}{"$start:$end"}=0;
	}
}

close BED;

while(<OUT_BED>){
    my ($chr,$start,$end,$info) = split("\t");
    my $bin1 = int($start/10000);
	my $bin2 = int($end/10000);
	foreach my $k($bin1..$bin2){
		$out_hash{$chr}{$k}{"$start:$end"}=0;
	}
}
close OUT_BED;

my $totalReads=0;
my $captured=0;
my $less_150=0;
my $less_500=0;
my $more_500=0;
while(<READ_pos>){
    chomp;
    my ($chr,$start,$end) = split("\t");
    my $bin = int($start/10000);
    $totalReads++;
########################################################
#Exome:	-----------@********************O---------------
#Reads:   @*******O								    	#no_captured
#Reads:    @*******O									#no_captured
#Reads:     @*******O								    #captured
#Reads:                    @*******O					#captured
#Reads:                                @*******O		#captured
#Reads:                                 @*******O		#no_captured
#Reads:  								 @*******O		#no_captured
########################################################
	my $flag=0;
    foreach my $key (keys %{$in_hash{$chr}{$bin}}){
        my ($exon_start, $exon_end) = split(":",$key);
		my $actual_exon_start = $exon_start + 1;
		my $actual_exon_end = $exon_end + 1;
        unless ($start >= $actual_exon_end or $end <= $actual_exon_start){
            $captured++;
			$flag=1;
            last;
        }
    }
########################################################
#Out_Exome:	-----------@********************O---------------
#Reads:      @*******O								  	#no_captured
#Reads:    	        @*******O							#no_captured
#Reads:                @*******O					    #captured
#Reads:                    @*******O					#captured
#Reads:                             @*******O			#captured
#Reads:                                 @*******O		#no_captured
#Reads:  				                    @*******O	#no_captured
########################################################
	if($flag==0){
		foreach my $key (keys %{$out_hash{$chr}{$bin}}){
			my ($region_start, $region_end) = split(":",$key);
			my $actual_region_start = $region_start + 1;
			my $actual_region_end = $region_end + 1;
			if	($start >= $actual_region_start and $end <= $actual_region_end){
				my $left = $start - $actual_region_start;
				my $right = $actual_region_end - $end;
				if($left<150 || $right<150){
					$less_150++;
				}elsif($left>500 && $right>500){
					$more_500++;
				}else{
					$less_500++;
				}
				$flag=2;
            	last;
			}
        }
	}
}
close READ_pos;


printf ("MappingReads\tCapturedReads\tless 150bp\tless 500bp\tmore 500bp\n%d\t%d\t%d\t%d\t%d\n",$totalReads,$captured,$less_150,$less_500,$more_500);
#system("date");
