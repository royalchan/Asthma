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
open(OUT, "> AA_RNA_Hits_Reordered_QN.txt") || die;

# First index the hits.
my %hit;
while ($line = <IN1>) {
	chomp $line;
	$hit{$line} = 1;
}

# Then stream the AA RNA-Seq raw data file.
$line = <IN2>; # Get header.
print OUT $line;
while ($line = <IN2>) {
	@linearray = split("\t", $line);
	my @linearray0 = split("_", $linearray[0]);
	my $key = $linearray0[0]."_".$linearray0[1];
	my $l = scalar @linearray0;
	next if (not defined $hit{$key} || 0);
	my $chr = join("_", @linearray0[2..($l-3)]);
	my $coordstart = $linearray0[$l - 2];
	my $coordend = $linearray0[$l - 1];
	$linearray[0] = $key."^".$chr."^".$coordstart."^".$coordend;
	my $reline = join("\t", @linearray);
	print OUT $reline;
}

close IN1;
close IN2;
close OUT;

exit;
