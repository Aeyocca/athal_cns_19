#! /usr/bin/perl -w

# bed_density.pl
# calculate density of elements in bed_b around the elements in bed_a
# Alan E. Yocca
# 06-04-19

use strict;
use Getopt::Long;

my $usage = "\n$0\n" .
	"\tAssumes input is sorted, does not check to see it is.. \n" .
	"\t--bed_a\n" .
		"\t\t<will output single value for EACH of these elements>\n" .
	"\t--genomecov\n" .
		"\t\t<genome coverage file (bed_b) of what you are calculating the density of>\n" .
	"\t--window\n" .
		"\t\t<number (in bp) to check around elements (upstream and downstream) in bed_a>\n" .
		"\t\t<inefficient program, increasing this should increase run time>\n" .
		"\t\t<default: 1000 (ie 1kb up and 1kb down)>\n" .
	"\t--intrafeature\n" .
		"\t\t<calculate coverage of elements in bed_a instead of density around them>\n" .
	"DOES NOT KEEP TRACK OF STRAND INFORMATION (not important to me atm, might add later)\n" .
	"\t--output\n\n";
	
#example command to get genomecov output (bedgraph format):
#[yoccaala@dev-intel14 01_bed]$ while read genome; do base=$(basename ${genome} | sed "s/\.fasta//"); echo "Starting ${base}"; time bedtools genomecov -bg -i /mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/26_repeat/01_bed/${base}_tair_te_e4_sort.bed -g ${genome}.fai > /mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/26_repeat/01_bed/${base}_tair_te_e4_sort_gc_bg.bed; done < /mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/16_final_genomes/final_genome_all_no_rep_list.txt

my $bed_a;
my $genomecov;
my $window=1000;
my $intrafeature = '';
my $output;

GetOptions ( "bed_a=s" => \$bed_a,
  "genomecov=s" => \$genomecov,
  "window=s" => \$window,
  "intrafeature" => \$intrafeature,
  "output=s" => \$output
) or die "$usage\n";

if (! defined $bed_a || ! defined $output || ! defined $genomecov) {
	die "$usage";
}

#loading entire bed_b into hash? tough on memory, but wayy easier to write
#basically getting bedtools genomecov output here but doing it myself, must not be as efficient
#hmmm save time by just requesting bedtools genomecov as bed_b??
#can make it so doesn't need to calculate every site, do exist check, zero if not, simply add if se\o
#yupp, requesting bedtools genomecov

#load into hash, hmm bedgraph format??
my %genomecov;
my $array_index=0;
my $chrom='';
open (my $genomecov_fh, '<', $genomecov) || die "Cannot open the genomecov file: $genomecov\n\n";
while (my $line = <$genomecov_fh>) {
	chomp $line;
	my @line = split("\t",$line);
	#new chrom, reset array_index
	$array_index = 0 if ($line[0] ne $chrom);
	$chrom = $line[0];
	#hooray for multidimensional hash
	#hash, where key equal to array, array will be 2d:
	#0 - start ival, 1 - end ival, 2 - coverage
	${$genomecov{$line[0]}}[$array_index] = [$line[1],$line[2],$line[3]];
	$array_index+=1;
}
close $genomecov_fh;

#load bed A into hash
my %bed_a;
$array_index=0;
open (my $bed_a_fh, '<', $bed_a) || die "Cannot open the bed_a file: $bed_a\n\n";
while (my $line = <$bed_a_fh>) {
	chomp $line;
	my @line = split("\t",$line);
	$array_index = 0 if ($line[0] ne $chrom);
	$chrom = $line[0];
	#chromosome, start, stop
	${$bed_a{$line[0]}}[$array_index] = [$line[1],$line[2]];
	$array_index+=1;
}
close $bed_a_fh;

#collapsing, fun part
my %output; #window_density, intragenic density, hash of arrays
#dont really need feature, hmm yes we do because need to ensure incrementing same feature
#name each feature chrom_start_stop?, equal to array?
my @splice;
my $last_index;

foreach my $chrom (sort keys %genomecov) {
	#foreach chromosome,
	#@{$hash{$string}}
	next if (! defined $bed_a{$chrom});
	#sort bed_a matching this chromosome into multidimensional array, well, two dimensional, start, stop
	my @sort_bed_a=();
	#for each feature on the specified chromosome
	#my $size = scalar @{$bed_a{$chrom}};
	#die "$size\n";
	for (my $interval_index = 0; $interval_index < scalar @{$bed_a{$chrom}}; $interval_index++) {
		my @tmp = ($bed_a{$chrom}[$interval_index][0],$bed_a{$chrom}[$interval_index][1]);
		#$bed_a{$chrom}[$interval_index][0]
		push @sort_bed_a, \@tmp;
	}
	#my $size = scalar @sort_bed_a;
	#die "$size\n";
	for (my $interval_index = 0; $interval_index < scalar @{$genomecov{$chrom}}; $interval_index++) {
		COORD:	for (my $coord = ${$genomecov{$chrom}}[$interval_index][0]; 
				$coord <= ${$genomecov{$chrom}}[$interval_index][1];
				$coord++) {
			#foreach coordinate in that interval
			#check bed_a hash and add there if appropriate
			#want to check the first few elements of bed_a, but dont want to sort every coordinate...
			#think will sort into multidimensional array..
		
			#splice things out, splice from end, thanks SO
			for ( sort { $b <=> $a } @splice ) {
	        	splice @sort_bed_a, $_, 1;
		    }
			@splice=();
			INDEX: for my $index (0 .. ((scalar @sort_bed_a) - 1)) {
				#name feature for output purposes
				#my $size = scalar( @{$sort_bed_a[$index]});
				#die "size: $size\ti-0: $sort_bed_a[$index][0]\n";
				my $feature = $chrom . "_" . $sort_bed_a[$index][0] . "_" . $sort_bed_a[$index][1];
				if (! defined $output{$feature}) {
					#first time encountering this bed_a feature, set window and intrafeature to zero
					#print "new def index\n";
					$output{$feature}[0] = 0;
					$output{$feature}[1] = 0;
				}
				#if not close enough, next coord
				if ($coord < ($sort_bed_a[$index][0] - $window)) {
					#capture rest of bed_a features below
					$last_index = $index;
					next COORD 
				}		
				
				#delete $index (bed_a feature) if past it, then next index
				#ugh that will skip one.., ahh splice after breaking, nice
				#just mark for splicing
				if ($coord > ($sort_bed_a[$index][1] + $window)) {
					push @splice, $index;
					next INDEX;
				}
				#check if intergenic or in window, above will break loop if outside window
				if ($coord < $sort_bed_a[$index][0]) {
					$output{$feature}[0]+=1 ;
				}
				elsif ($coord > $sort_bed_a[$index][1]){
					$output{$feature}[0]+=1;
				}
				elsif ($coord >= $sort_bed_a[$index][0] && $coord <= $sort_bed_a[$index][1]) {
					#must be intrafeature.. alright check
					$output{$feature}[1]+=1;
				}
				else {
					die "Something went wrong:\n" . "coord: $coord\n" .
					"start: $sort_bed_a[$index][0]\n" . "stop: $sort_bed_a[$index][1]\n";
				}
			}
		}
	}
	#finish splicing, loop through rest of bed_a features, setting things to zero
	#ooo, don't need to splice if starting at $last_index
#	for ( sort { $b <=> $a } @splice ) {
#	    splice @sort_bed_a, $_, 1;
#	}
#	@splice=();
	#Start at $last_index+1 to get features not yet set to zero, eh, set the feature at last_index = 0 already anyway
	#my $size = scalar @sort_bed_a;
	#die "last_index: $last_index\nsize: $size\n";
	for my $index ($last_index .. ((scalar @sort_bed_a) - 1)) {
		my $feature = $chrom . "_" . $sort_bed_a[$index][0] . "_" . $sort_bed_a[$index][1];
		$output{$feature}[0] = 0;
		$output{$feature}[1] = 0;
		#print "new def outside\n";
	}
}

open (my $out_fh, '>', $output) || die "Cannot open the output file: $output\n\n";
foreach my $key (keys %output) {
	#normalize feature size
	my @feature = split("_",$key);
	my $window_size = $window * 2;
	my $intra_size = $feature[2] - $feature[1] + 1;
	my $window_prop = $output{$key}[0] / $window_size;
	my $intra_prop = $output{$key}[1] / $intra_size;
	print $out_fh "$window_prop\tWindow\n";
	print $out_fh "$intra_prop\tIntrafeature\n";
}
close $out_fh;

exit;


