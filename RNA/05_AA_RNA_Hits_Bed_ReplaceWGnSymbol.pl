#!/usr/bin/perl

use strict;
use warnings;

####################################
# This script is to replace RefSeq Accession Numbers with Gene Symbols in the bed files.
####################################

use FileHandle;

my $line;
my @linearray;

open(FILES, "05_AA_RNA_Hits_Bed_ReplaceWGnSymbol_Files.txt") || die;
my @filenames;
my @outnames;
while ($line = <FILES>) {
	chomp $line;
	push (@filenames, $line);
	my $tmp = $line;
	$tmp =~ s/\.bed//;
	$tmp .= "_genesymbol.bed";
	push (@outnames, $tmp);
}
close FILES;

# Index the conversion table.
open(TABLE, "AA_RNA_Hits_AccNo_vs_GeneSymbol.txt") || die;
my %table;
while ($line = <TABLE>) {
	chomp $line;
	@linearray = split("\t", $line);	
	$table{$linearray[0]} = $linearray[1];
}
close TABLE;


my $fh = FileHandle -> new;
my $oh = FileHandle -> new;

for (my $i = 0; $i < (scalar @filenames); ++$i) {

	$fh -> open ("< $filenames[$i]");
	$oh -> open ("> $outnames[$i]");
	
	# Output header.
	$line = <$fh>;
	print $oh $line;
	
	# convert and print.
	while ($line = <$fh>) {
	
		chomp $line;	
		@linearray = split("\t", $line);
		
		my @linearray3 = split("\\^", $linearray[3]);
		if ($filenames[$i] =~ /Wider/) {
			$linearray[3] = $table{$linearray3[0]};
		} else {
			$linearray[3] = $table{$linearray3[0]}."\^".$linearray3[1];
		}
		
		my $reline = join ("\t", @linearray);
		print $oh "$reline\n";
		
	}
	
	$fh -> close;
	$oh -> close;
	
}

exit;
