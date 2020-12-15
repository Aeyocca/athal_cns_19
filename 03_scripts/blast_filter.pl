#! /usr/bin/perl -w

# blast_filter.pl
# take in blast file and make as many filters as you can for it
# original purpose: take blast between augustus annotations for accessions and TAIR10 proteome and assign an AT locus format geneID to the augustus annotations
# do accomplish original purpose, we need best ref hit for every query
# Alan E. Yocca
# 09-10-18

use strict;
use warnings;
use Getopt::Long;


my $usage = "\n$0\n" .
			"\t--blast <blast input m8 format>\n" .
			"\t--number_hits <keep up to the top x number of blast hits based on specified criteria>\n" .
				"\t\t<default: 5>\n" .
				"\t\t<set to insanely high value to keep all>\n" .
			"\t--queryid <sort by query id>\n" .
				"\t\t<2 POSSIBLE OPTIONS:>\n" .
					"\t\t\t<--queryid character (or really any word starting with the letter c)>\n" .
						"\t\t\t\t<sort as a character>\n" .
					"\t\t\t<--queryid numeric (or any word starting with the letter n)>\n" .
						"\t\t\t\t<sort as numeric>\n" .
			"\t--refid <sort by ref id (see queryid for options)>\n" .
			"\t--bit_score <bitscore value above which to keep>\n" .
				"\t\t<default: 1>\n" .
			"\t--ref_start <start of ref hits to cut>\n" .
				"\t\t<if multiple hits have the same start, will sort by longest to shortest hit>\n" .
			"\t--query_start <start of query hits to cut>\n" .
				"\t\t<if multiple hits have the same start, will sort by longest to shortest hit>\n" .
			"\t--id <id percentage to filter>\n" . 
			"\t--align_length <alignment length to filter out below>\n" .
			"\t--align_percent\n" .
				"\t\t<make sure blast hits greater than X% length of query sequence>\n" .
				"\t\t<EXPRESS AS DECIMAL>\n" .
				"\t\t<IF SPECIFIED, required to specify query fasta as well>\n" .
			"\t--query_fasta\n" .
				"\t\t<if align_percent specified>\n" .
			"\t--e_value <e-value to filter out below>\n" .
				"\t\t<default 1e-4>\n" .
			"\t--filter_order <comma separated list ordered of things to filter out, match command line flags exactly>\n" .
				"\t\t<optional, if you want, if not defined, will just take the order of things on command line>\n" .
				"\t\t<default: bit_score, then alignment length>\n" .
			"\t--no_self\n" .
				"\t\t<only report non-self hits>\n" .
			"\t-o <output filtered blast> \n\n";

#filter brainstorm
#bit score
#length ref
#length query
#id %
#alignment length
#mismatches (not really, why would want to do this?)
#gap openings (same as mismatches)
#e-value

#also, how to programmatically decide order, also need thresholds for these
#easiest would be order on cmd line right? could also make another flag for it
#call GetOptions a bunch of times, yupp
#no, can do $1, $2, etc though

my $refid;
my $align_percent;
my $query_fasta;
my $force='';
my $queryid;
my $blast;
my $output;
my $filter_order;
my $number_hits=5;
my $bit_score=1;
my $ref_start;
my $query_start;
my $id;
my $align_length;
my $e_value=0.0001;
my $no_self = '';

GetOptions ( "blast=s" => \$blast,
  "force" => \$force,
  "refid=s" => \$refid,
  "queryid=s" => \$queryid,
  "o=s" => \$output,
  "filter_order=s" => \$filter_order,
  "number_hits=i" => \$number_hits,
  "bit_score=i" => \$bit_score,
  "ref_start=i" => \$ref_start,
  "query_start=i" => \$query_start,
  "id=i" => \$id,
  "align_length=i" => \$align_length,
  "align_percent=s" => \$align_percent,
  "query_fasta=s" => \$query_fasta,
  "no_self" => \$no_self,
  "e_value=s" => \$e_value
) or die "Not load opts\n$usage\n";

if ( (!(defined $blast)) || (!(defined $output)) ) {
  print "$usage\n" . "No blast or output defined\n";
  exit;
}

if ( (defined $refid) ) {
	if ( $refid !~ /^[nc]/ ) {
		die "$usage\n" . "refid defined but not in bounds\n";
	}
}

if ( (defined $queryid) ) {
	if ( $queryid !~ /^[nc]/ ) {
		die "$usage\n" . "queryid defined but not in bounds\n";
	}
}

if (defined $align_percent && ! defined $query_fasta) {
	print "$usage\n";
	die "If align_percent given, need to specify query_fasta\n";
}

if (-e $output && (!(defined $force))) {
print "File: $output or exist, is it okay to overwrite it?\n"; 
my $answer = <STDIN>;
	if ($answer =~ /^y(?:es)?$/i) {
		print "Excellent!\n";
	}
	else {
		die "fine, I will not overwrite your files, but I will also not run this script\n";
	}
}

if (!(defined $filter_order)) {
	#initialize with comma so can add to it in the loop
	my $temp_filter_order = ",";
	for (my $i=0; $i<@ARGV; $i++) {
		next if ($ARGV[$i] !~ /^--/);
		next if ($ARGV[$i] eq "blast" || $ARGV[$i] eq "number_hits");
		$temp_filter_order = $temp_filter_order . $ARGV[$i] . ",";
	}
	$filter_order = $temp_filter_order;
}

#hopefully it will work
#test
#print "$filter_order\n";

open (my $blast_fh, '<', $blast) || die "Cannot open the blast file $blast\n\n";
open (my $out_fh, '>', $output) || die "Cannot open output fasta: $output\n\n";

my %query_fasta_length;
my $header;
if (defined $query_fasta) {
	print "collecting query sequence lengths\n\n";
	open (my $qf_fh, '<', $query_fasta) || die "Cannot open the query fasta file $query_fasta\n\n";
	while (my $line = <$qf_fh>) {
		chomp $line;
		if ($line =~ /^>/) {
			$line =~ s/^> *//g;
			#separate from $header variable
			my @header = split(" ", $line);
			#takes everything before the first space........
			#what if space right after carrot...
			#fixed above
			$header = $header[0];
		}
		else {
			$query_fasta_length{$header} += length($line);
		}
	}
	close $qf_fh;
}

#blast m8 format with array index:
#0: query	
#1: subject	
#2: %id	
#3: alignment length
#4: mismatches	
#5: gap openings	
#6: query start
#7: query end
#8: subject start
#9: subject end	
#10: E value
#11: bit score

#actually, should just load cns into hash, then loop through blast


#load blast into 2d array to sort, 
#this might get a little dicey if big blast file, 
#I'm working with arabidopsis and on a cluster so I am not worried, 
#not taking the time to do otherwise, 
#I mean, you have to load it all in to sort, right?


my @blast;
my $skipped=0;
LOAD_BLAST: while (my $line = <$blast_fh>) {
	chomp $line;
	next if ($line =~ /^#/);
	my @info = split("\t",$line);
	next if ($no_self && $info[0] eq $info[1]);
	next if ($info[10] > $e_value);
	push(@blast, \@info);
}

my $blast_lines = scalar @blast;
print "Loaded blast lines: $blast_lines\n";

my @filter_order = split(",",$filter_order);

my $filter_number = 0;
#define sort string
my $sort_code= " ";
#for the queryid / refid sorting:
my $op;

for (my $i=0; $i<@filter_order; $i++) {
	#not even sure this will work, try test
	#yea this isn't going to workout
	#make strings to put.. but how can interpret variable in string.... maybe different subroutines depending on number of fields then separate strings for the number of fields, then would still be six different loops, still a little different if starts is one.. we can try
	#can we just sort by all columns, but the b a thing switching on start / stop.. :/
	next if (!defined $filter_order[$i]);
	if ($filter_order[$i] eq "bit_score") {
		$sort_code = $sort_code . '\$b->[11] <=> \$a->[11]' . ' || ';
	}
	elsif ($filter_order[$i] eq "ref_start") {
		$sort_code = $sort_code . '\$a->[8] <=> \$b->[8]' . ' || ';
		$sort_code = $sort_code . '\$b->[9] <=> \$a->[9]' . ' || ';
	}
	elsif ($filter_order[$i] eq "query_start") {
		$sort_code = $sort_code . '\$a->[6] <=> \$b->[6]' . ' || ';
		$sort_code = $sort_code . '\$b->[7] <=> \$a->[7]' . ' || ';
	}
	elsif ($filter_order[$i] eq "id") {
		$sort_code = $sort_code . '\$b->[2] <=> \$a->[2]' . ' || ';
	}
	elsif ($filter_order[$i] eq "align_length") {
		$sort_code = $sort_code . '\$b->[3] <=> \$a->[3]' . ' || ';
	}
	elsif ($filter_order[$i] eq "e_value") {
		$sort_code = $sort_code . '\$a->[10] <=> \$b->[10]' . ' || ';
	}
	elsif ($filter_order[$i] eq "refid") {
		if ($refid =~ /^n/) {
			$op = '<=>';
		}
		else {
			$op = 'cmp';
		}
		$sort_code = $sort_code . '\$a->[1] $op \$b->[1]' . ' || ';
	}
	elsif ($filter_order[$i] eq "queryid") {
		if ($queryid =~ /^n/) {
			$op = '<=>';
		}
		else {
			$op = 'cmp';
		}
		$sort_code = $sort_code . '\$a->[0] $op \$b->[0]' . ' || ';
	}
	else {
		die "Make sure the filter_order matches one of the optional flags\n$usage";
	}
	$filter_number = $filter_number + 1;
}

#defaults
my @sorted_blast;
if ($filter_number == 0) {
	#default
	print "using default sort: bit score then alignment length\n";
	$sort_code = $sort_code . '\$a->[11] <=> \$b->[11]' . ' || ';
	$sort_code = $sort_code . '\$a->[3] <=> \$b->[3]' . ' || ';
	@sorted_blast = sort {
		$b->[11] <=> $a->[11] ||
		$b->[3] <=> $a->[3]
	} @blast;
}
else {
	#remove last ' || ' added
	$sort_code =~ s/\s\|\|\s$//;
	@sorted_blast = sort {custom_sort()} @blast;
	sub custom_sort {
		eval qq{$sort_code};
	}
}

my $sort_blast_lines = scalar @sorted_blast;
print "Loaded sorted blast lines: $sort_blast_lines\n";

#test
#print "$sort_code\n";
#die;

#sort blast 2d array, hmm, should group by query hit right? like queryID
#ahh I think I want to sort each is what I wants, best way to do that? first sort by query ID, yupp, well thats for my purposes, I'll make sure to specify on cmd line



#hmm, just print out now right? ahh make hash to trim if number of hits falls
#not sure why you would want to filter out after x amount of specific reference hits, if you really need that, either rerun the blast search, or simply switch the ref / query columns in your blast file

#hash with queryid => number of times seen
my %seen_it;
my $dropped_length=0;

PRINT_BLAST: for (my $i=0; $i<@sorted_blast; $i++) {
	if ($seen_it{$sorted_blast[$i][0]}) {
		if ($seen_it{$sorted_blast[$i][0]} >= $number_hits) {
			next PRINT_BLAST;
		}
		else {
			$seen_it{$sorted_blast[$i][0]} = $seen_it{$sorted_blast[$i][0]} + 1;
		}
	}
	else {
		#set to 1
		$seen_it{$sorted_blast[$i][0]} = 1;
	}
	if (defined $align_percent) {
		#compute alignment length
		my $percent = $sorted_blast[$i][3] / $query_fasta_length{$sorted_blast[$i][0]};
		if ($percent < $align_percent) {
			$dropped_length +=1;
			next PRINT_BLAST;
		}
	}
	else {
		#nothing
	}
	for (my $j=0; $j < @{$sorted_blast[$i]}; $j++) {
		print $out_fh "$sorted_blast[$i][$j]\t";
	}
	print $out_fh "\n";
}

print "Dropped because of align_percent:\t$dropped_length\n" if (defined $align_percent);

close $out_fh;
close $blast_fh;

exit;


