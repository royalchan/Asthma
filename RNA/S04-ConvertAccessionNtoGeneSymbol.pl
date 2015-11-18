#!/usr/bin/perl

use strict;
use warnings;

####################################
#	This script is to swap the accession numbers with gene symbols in the AA37TxExp_Reordered_QN.txt file.
####################################

open(IN, "AA37TxExp_Reordered_QN.txt") || die;
open(OUT, "> AA37TxExp_Reordered_QN_GeneSymbol.txt") || die;
open(CONVERSIONTABLE, "AA37TxExp_RefSeqIDtoGeneSymbol_biomaRt.txt") || die;
open(CONVERSIONTABLE2, "AA37TxExp_RefSeqIDtoGeneSymbol_DAVID.txt") || die;
open(CONVERSIONTABLE3, "AA37TxExp_RefSeqIDtoGeneSymbol_SOURCE.txt") || die;

my $line;
my @linearray;

$line = <CONVERSIONTABLE>;
my %ct;
while ($line = <CONVERSIONTABLE>) {
	chomp $line;
	@linearray = split("\t", $line);
	$ct{$linearray[0]} = $linearray[1];
}
while ($line = <CONVERSIONTABLE2>) {
	chomp $line;
	@linearray = split("\t", $line);
	$ct{$linearray[0]} = $linearray[1];
}
while ($line = <CONVERSIONTABLE3>) {
	chomp $line;
	@linearray = split("\t", $line);
	$ct{$linearray[0]} = $linearray[1];
}

$line = <IN>;
print OUT $line;
while ($line = <IN>) {
	@linearray = split("\t", $line);
	my @linearray0 = split("_", $linearray[0]);
	my $key = $linearray0[0]."_".$linearray0[1];
	if ((defined $ct{$key} && $ct{$key} ne "" && $ct{$key} ne "Data not found") || 0) {
		$linearray[0] = $ct{$key};
	} else {
		$linearray[0] = $key;	
	}
	my $reline = join("\t", @linearray);
	print OUT $reline;
}

close IN;
close OUT;
close CONVERSIONTABLE;
close CONVERSIONTABLE2;

exit;
