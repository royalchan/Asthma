#!/usr/bin/perl

use strict;
use warnings;

####################################
# This script is to generate 2kb and 5 kb bed files for the AA RNA-Seq hits.
####################################

my $line;
my @linearray;

open(IN1, "AA_RNA_Hits.txt") || die;
open(IN2, "AA40TxExp_Reordered_QN.txt") || die;
open(OUT1, "> AA_RNA_Hits_2kbWider.bed") || die;
open(OUT2, "> AA_RNA_Hits_5kbWider.bed") || die;

# First index the hits.
my %hit;
while ($line = <IN1>) {
	chomp $line;
	$hit{$line} = 1;
}

# Then stream the AA RNA-Seq raw data file.
$line = <IN2>; # Get rid of header.
print OUT1 "track name=\"aa_rna_hits_2kb\" description=\"AA RNA Hit Coordinates +- 2 kb\" visibility=dense\n";
print OUT2 "track name=\"aa_rna_hits_5kb\" description=\"AA RNA Hit Coordinates +- 5 kb\" visibility=dense\n";
while ($line = <IN2>) {
	@linearray = split("\t", $line);
	my @linearray0 = split("_", $linearray[0]);
	my $key = $linearray0[0]."_".$linearray0[1];
	my $l = scalar @linearray0;
	next if (not defined $hit{$key} || 0);
	my $chr = join("_", @linearray0[2..($l-3)]);
	my $coordstart1 = $linearray0[$l - 2] - 2000;
	my $coordstart2 = $linearray0[$l - 2] - 5000;
	my $coordend1 = $linearray0[$l - 1] + 2000;
	my $coordend2 = $linearray0[$l - 1] + 5000;
	print OUT1 "$chr\t$coordstart1\t$coordend1\t$key\n";
	print OUT2 "$chr\t$coordstart2\t$coordend2\t$key\n";
}

close IN1;
close IN2;
close OUT1;
close OUT2;

exit;
