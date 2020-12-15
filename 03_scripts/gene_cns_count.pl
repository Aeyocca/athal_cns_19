#! /usr/bin/perl -w

# gene_cns_count.pl
# cns_gene_assignment only looks within a certain window
# also it will double count CNS
# This is relatively dummy script.. yea it is because need specifically:
#	- output from bedtools closest assigning each CNS their proximate gene
# output will be the same as cns_gene_assignment.pl though so can just spit back into R code
# ripped most code from closest_bed_normalization.pl
# hmmmm probs easier to just rerun bedtools closest on mx_all annotations instead of having to load in multiple files? wait, just concatenate the coll and moved
# how to annotate though? wait doesn't matter for this, not worth investigating
# Alan E. Yocca
# 03-27-19

use strict;
use Getopt::Long;

my $usage = "\n$0\n" .
				"\t--closest_bed\n" .
				"\t--gene_list\n" .
					"\t\t<optional, if want to count zeros>\n" .
				"\t--fasta\n" .
					"\t\t<instead of gene_list, feed fasta file because I am lazy " .
					"and don't want to make a separate gene_list file>\n" .
				"\t--strip_transcript\n" .
					"\t\t<fasta file has transcript ids, take those off because they are not in bed file>\n" .
				"\t--gene_column\n" .
					"\t\t<array index of column containing genes that we count in --closest_bed file>\n" .
					"\t\t<default: 9>\n" .
				"\t--band_aid\n" .
					"\t\t<multiple matching gene names.. ugh, take positional information into account>\n" .
				"\t--posv_file\n" .
					"\t\t<specify PosV CNS file, will assign every collinear CNS" .
					" (those not here but found in bed file) to gene based on name not location>\n" .
#				"\t--trans_file\n" .
#						"\t\t<needed if --posv_file specified>\n" .
				"\t--no_meta_info\n" .
				"\t--output\n\n";
				
my $closest_bed;
my $gene_list;
my $fasta;
my $output;
my $strip_transcript = '';
my $gene_column = 9;
my $band_aid = '';
my $posv_file;
my $trans_file;
my $no_meta_info='';

GetOptions (
	"closest_bed=s" => \$closest_bed,
	"gene_list=s" => \$gene_list,
	"fasta=s" => \$fasta,
	"strip_transcript" => \$strip_transcript,
	"gene_column=s" => \$gene_column,
	"band_aid" => \$band_aid,
	"posv_file=s" => \$posv_file,
	"trans_file=s" => \$trans_file,
	"no_meta_info" => \$no_meta_info,
	"output=s" => \$output
) or die "$usage\n";

if (! defined $closest_bed || ! defined $output ) {
	die "Bed or output not defined\n$usage\n";
}
if (defined $gene_list && defined $fasta) {
	die "Specify one of --gene_list / --fasta\n$usage";
}

die "Band aid function wont work, better hope you have unique gene_ids\n" if ($band_aid);

open (my $cb_fh, '<', $closest_bed) || die "Cannot open the closest bed file: $closest_bed\n\n";
#0: Chr1
#1:    4931
#2:    4957
#3:    1||30332676||30332702||18879_AT1G80700||+||CNS||1||16157
#4:        500
#5:     +
#6:       Chr1
#7:    4826
#8:    5488
#9:    AT1G80980
#10:       500
#11:     +
#12:       0

#gene hash just going to straight up count should be quick
my %gene_cns_count;

if (defined $gene_list) {
	open (my $gl_fh, '<', $gene_list) || die "Cannot open the gene_list file: $gene_list\n\n";
	while (my $line = <$gl_fh>) {
		chomp $line;
		$gene_cns_count{$line} = 0;
	}
	close $gl_fh;
}
if (defined $fasta) {
	open (my $fasta_fh, '<', $fasta) || die "Cannot open the fasta file: $fasta\n\n";
	while (my $line = <$fasta_fh>) {
		chomp $line;
		next if ($line !~ /^>/);
		( my $transcript = $line ) =~ s/^> *//g;
		my @transcript = split(" ",$transcript);
		my $gene;
		if ($strip_transcript) {
			( $gene =  $transcript[0] ) =~ s/\.[0-9]*//g;
		}
		else {
			$gene = $transcript[0];
		}
		#die "$gene\n";
		$gene_cns_count{$gene} = 0;
	}
	close $fasta_fh;
}

my %trans;
my %posv_cns;
if (defined $trans_file && defined $posv_file) {
#if (defined $posv_file) {
	#have a version where don't need to translate
	open (my $trans_file_fh, '<', $trans_file) || die "Cannot open the trans_file : $trans_file\n\n";
	while ( my $line = <$trans_file_fh>) {
		chomp $line;
		my @line = split("\t",$line);
		$line[1] =~ s/\..*//g;
		#die "$line[1]\n";
		$trans{$line[1]} = $line[0];
	}
	close $trans_file_fh;
	open (my $posv_file_fh, '<', $posv_file) || die "Cannot open the posv_file : $posv_file\n\n";
	while ( my $line = <$posv_file_fh>) {
		chomp $line;
		$posv_cns{$line} = 1;
	}
	close $posv_file_fh;
}



#go through bed
my %coll_cns_count;
my $bed_line_count = 0;
my $coll_count = 0;
BED: while (my $line = <$cb_fh>) {
	chomp $line;
	my @line = split("\t",$line);
	#non-zero check so gene_list / fasta can be optional
	my $gene = $line[$gene_column];	
	if ($band_aid) {
		#gene is now: start_stop_gene
		$gene = $line[7] . "@" . $line[8] . "@" . $line[9];
		#delete $gene_cns_count{$line[9]} if (defined $gene_cns_count{$line[9]});
	}
	
	##############
	#count towards the gene assignment if not a posV CNS
	#die "$line[3]\n";
	my @meta_info = split(/\|\|/,$line[3]) if (! $no_meta_info);
	my $cns = ($no_meta_info) ? $line[3] : $meta_info[3];
	if ( ! defined $posv_cns{$cns}  && defined $posv_file) {
		#assign based on CNS name, not $gene
		my @gene = split("_",$cns);
		#$trans{$gene[1]}
		#die "$gene[1]\n";
		if (! defined $trans{$gene[1]}) {
			#the gene this CNS was associated with is not in trans file, so not collinear anymore, assign to proximate gene??????????????
			#yea, think thats how it would function if CNS loses gene association right??
			#skip out of this if loop
			if (defined $gene_cns_count{$gene} && $gene_cns_count{$gene} != 0) {
				$gene_cns_count{$gene}+=1;
				#die "$gene\n" if ($gene =~ /AT/);
			}
			else {
				$gene_cns_count{$gene}=1;
				#die "$gene\n" if ($gene =~ /AT/);
			}
			$bed_line_count+=1;
			next BED;
		}
		#die "defined: $trans{$gene[1]}\n";
		if (defined $gene_cns_count{$trans{$gene[1]}} && $gene_cns_count{$trans{$gene[1]}} != 0) {
			$gene_cns_count{$trans{$gene[1]}}+=1;
			#die "$trans{$gene[1]}\n" if ($trans{$gene[1]} =~ /AT/);
		}
		else {
			$gene_cns_count{$trans{$gene[1]}}=1;
			#die "$trans{$gene[1]}\n" if ($trans{$gene[1]} =~ /AT/);
		}
		$bed_line_count+=1;
		$coll_count+=1;
		next BED;
	}
	##############
	
	
	if (defined $gene_cns_count{$gene} && $gene_cns_count{$gene} != 0) {
		$gene_cns_count{$gene}+=1;
		#die "$gene\n" if ($gene =~ /AT/);
	}
	else {
		$gene_cns_count{$gene}=1;
		#die "$gene\n" if ($gene =~ /AT/);
	}
	$bed_line_count+=1;
}
#print "coll count: $coll_count\n";
#print "Bed count: $bed_line_count\n";
#need one more bandaid fix, add collinear counts to each duplicate

#print out
#output example:
#Gene    CNS_Count       Type
#augustus_masked-Chr1-abinit-gene-0.0-mRNA-1     0       Conserved
#augustus_masked-Chr1-abinit-gene-0.1-mRNA-1     2       Conserved
#augustus_masked-Chr1-abinit-gene-0.9-mRNA-1     1       Conserved
#drop type...

open (my $out_fh, '>', $output) || die "Cannot open the output file: $output\n\n";

print $out_fh "Gene\tCNS_Count\n"; 
#my $size = keys %gene_cns_count;
#print "size gene_cns_count: $size\n";

my $print_count = 0;
PRINT: foreach my $key (keys %gene_cns_count) {
	my $gene = $key;
	if ($band_aid) {
		if ($gene =~ /@/) {
			#posv only counted for this gene, we need 
			#to add the collinear counts
			$gene =~ s/.*@//g;
			#die "$gene\n";
			$gene_cns_count{$key} += $gene_cns_count{$gene} if (defined $gene_cns_count{$gene});
			print $out_fh "$gene\t$gene_cns_count{$key}\n";
			$print_count+=1;
			next PRINT;
		}
		#else print as usual
	}
	#die "got out\n";
	$print_count+=1;
	print $out_fh "$gene\t$gene_cns_count{$key}\n";
}

#print "p coutn: $print_count\n";
	
	
#		print $out_fh "$gene\t$gene_cns_count{$key}\n" if ($gene !~ /_/);
##		if ($gene !~ /_/) {
##			die "grep worked: $gene\n";
##		}
#		next PRINT if ($gene !~ /_/);
#		$gene =~ s/.*_//g;
#		#that should fix it, 
#		$gene_cns_count{$key} += $gene_cns_count{$gene} if (defined $gene_cns_count{$gene});
#		if (defined $gene_cns_count{$gene}) {
#			die "Added\n";
#		}
#	}
#	print $out_fh "$gene\t$gene_cns_count{$key}\n";
#}

close $out_fh;














