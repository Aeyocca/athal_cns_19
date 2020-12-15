#! /usr/bin/perl -w

# prox_feature_distance.pl
# take in output of bedtools closest, for each line:
# calculate distance of feature A to feature B
# added function to take in list? or just wanna make separate file for each ecotype?
# separate file, then can read them all in in R
# Alan E. Yocca
# 05-02-19

use strict;
use Getopt::Long;

my $usage = "\n$0\n" .
	"\t--input\n" .
	"\t--header\n" .
		"\t\t<no arg, if specified will add header>\n" .
	"\t--summary_output\n" .
		"\t\t<kind of a dummy flag>\n" .
		"\t\t<specify file to append avg / stdev distance to prox feature>\n" .
	"\t--quiet\n" .
	"\t--output\n\n";

my $input;
my $output;
my $header = '';
my $summary_output;
my $quiet = '';

GetOptions ( "input=s" => \$input,
  "header" => \$header,
  "summary_output=s" => \$summary_output,
  "quiet" => \$quiet,
  "output=s" => \$output
) or die "$usage\n";

if (! defined $input ) {
	die "$usage";
}

open (my $input_fh, '<', $input) || die "Cannot open the input file: $input\n\n";
open (my $out_fh, '>', $output) || die "Cannot open the output file: $output\n\n" if (defined $output);

print $out_fh "Feature_A\t" ."Feature_B\t" . "Distance\t" . 
					"Feature_A_Strand\t" . "Feature_B_Strand\n" if (defined $out_fh && $header);

my @dist_vect;

while (my $line = <$input_fh>) {
	chomp $line;
	my @line = split("\t",$line);
	#example line:
	#0 - Chr1    
	#1 - 11074   
	#2 - 11129   
	#3 - 1||11153||11098||37061_AT1G01030||-||CNS||1||14_1       
	#4 - 500     
	#5 - -       
	#6 - Chr1    
	#7 - 11840   
	#8 - 12916   
	#9 - AT1G01030    
	#10 - 500     
	#11 - -       
	#12 - 712
	
	#see if upstream, downstream, or inter. if inter just putting zero
	#five scenarios to check, hmm can make all intergenic single so 3 scenarios: fudge need to think about strand..
	#do up/down based on gene strand
	my $flip = ($line[11] eq "-" ? "true" : "false");
	
	#run through plus strand scenario then see how to most efficiently handle minus strand
	my $ltss_dist = $line[1] - $line[7];
	my $rtss_dist = $line[2] - $line[7];
	my $ltes_dist = $line[1] - $line[8];
	my $rtes_dist = $line[2] - $line[8];
	
	my $dist;
	if ( ($ltss_dist * $rtss_dist < 0 ) || ($ltes_dist * $rtes_dist < 0 ) ) {
		#intergenic, straddles gene boundary
		$dist = 0;
	}
	elsif ( ($ltss_dist < 0 ) & ( $rtss_dist < 0 ) ) {
		#left side
		#so subtracting from left side (tss) and taking smallest value
		$dist = (abs($ltss_dist) < abs($rtss_dist) ? $ltss_dist : $rtss_dist);
	}
	else {
		#right side
		#so subtracting from right side (tss) and taking smallest value
		$dist = (abs($ltes_dist) < abs($rtes_dist) ? $ltes_dist : $rtes_dist);
	}
	#flip sign of distance if gene is minus strand
	#before flipping, + strand upstream dist will be negative
	$dist *= -1 if ($flip eq "true");

	print $out_fh "$line[3]\t" . "$line[9]\t" . "$dist\t" . "$line[5]\t" . "$line[11]" . "\n" if (defined $out_fh);
	push @dist_vect, $dist;
}

close $input_fh;
close $out_fh if (defined $out_fh);

if (defined $summary_output) {
	print "Opening: $summary_output for appending\n" if (! $quiet);
	open (my $sout_fh, '>>', $summary_output) || die "Cannot open the summary_output file: $summary_output\n\n";
	my @abs_dist;
	foreach (@dist_vect) {
		my $tmp = sqrt($_ ** 2);
		push @abs_dist, $tmp;
	}
	my $avg = average(@abs_dist);
	my $std = stdev(@abs_dist);
	print $sout_fh "$avg\t$std\n";
}

exit;

sub average{
        my(@data) = @_;
        if (scalar @data == 0) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@data) {
                $total += $_;
        }
        my $average = $total / scalar(@data);
        return $average;
}
sub stdev{
        my(@data) = @_;
        if(scalar(@data) == 1){
                return 0;
        }
        my $average = average(@data);
        my $sqtotal = 0;
        foreach(@data) {
                $sqtotal += ( ($average-$_) ** 2 );
        }
        my $std = ($sqtotal / (scalar(@data)-1)) ** 0.5;
        return $std;
}

