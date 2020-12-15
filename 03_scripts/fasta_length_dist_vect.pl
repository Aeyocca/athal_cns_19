#! /usr/bin/perl -w

# fasta_length_dist_vect.pl
# take fasta file, output a comma separated vector (list??, its just one line) of the lengths of each fasta entry
# think of adding an ignore N switch at some point
# WILL NOT DO ANY CHECKS, simply get number of characters in fasta entry
# should make it work on wrapped fasta
# 06-28-18
# Alan E. Yocca
# 10-15-18
# add some more functionality, would like to take in bed file, hmm I guess that should be a separate script huh, 
# well I can just read bed file straight into R so don't even need a script
# alrighty, will add functionality to this later when / if I need it
# 12-12-18
# Changed to do first column header, second column (tab separated) length
# or just single column with length

use strict;
use Getopt::Long;

my $usage = "\n$0 -f <input fasta> --header <optional include header> -o <output> \n\n";

my $opt_f;
my $header_out;
my $opt_o;

GetOptions ( "f=s" => \$opt_f,
  "header" => \$header_out,
  "o=s" => \$opt_o
) or die "$usage\n";

if ( (!(defined $opt_f)) || (!(defined $opt_o)) ) {
  print "$usage";
  exit;
}

open (my $fasta_fh, '<', $opt_f) || die "Cannot open the fasta file: $opt_f\n\n";
open (my $out_fh, '>', $opt_o) || die "Cannot open the output file: $opt_o\n\n";

#my $unambig = 0;
#my $ns = 0;
my $fasta_entries = 0; #just a check to see if working properly

my %length;
my $header;

while (my $line = <$fasta_fh>) {
	chomp $line;
	if ($line =~ m/^>/) {
		#declare numeric so can add to it
		$length{$line} = 0;
		$header = $line;
	}
	else {
		$length{$header} += length($line);
	}
}

#print out
foreach my $fasta_header (keys %length) {
	if ($header_out) {
		print $out_fh "$fasta_header\t$length{$fasta_header}\n";
	}
	else {
		print $out_fh "$length{$fasta_header}\n";
	}
}

close $fasta_fh;
close $out_fh;

exit;
