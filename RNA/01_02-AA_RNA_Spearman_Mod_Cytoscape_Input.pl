#!/usr/bin/perl

use strict;
use warnings;

# This script is to generate assumed cytoscape correlation network input
#	from the AA37TxExp_Spearman_RhoSquared.txt and the AA37TxExp_Spearman_Rho.txt files.

open(IN, "AA_DA_RNA_Spearman_RhoSquared.txt") || die "Cannot open the input file!!\n\n";
open(OUT, "> AA_DA_RNA_Spearman_RhoSquared_CytoscapeInput.txt");

my $line;
my @linearray;

$line = <IN>; # Get Time Point Name Column.
chomp $line;
my @colnames = split("\t", $line);
shift(@colnames);
print OUT "AA_Sample_1\tAA_Sample_2\tSpearman_Rho_Squared\n";

my $i = 0;
while ($line = <IN>) {
	chomp $line;
	@linearray = split("\t", $line);
	my $rowname = shift(@linearray);
	for (my $j = $i + 1; $j < scalar @linearray; ++$j) {
		if ($linearray[$j] >= 0.16) {
			print OUT "$rowname\t$colnames[$j]\t$linearray[$j]\n";
		}
	}
	++$i;
}

close IN;
close OUT;

open(IN, "AA_DA_RNA_Spearman_Rho.txt") || die "Cannot open the input file!!\n\n";
open(OUT, "> AA_DA_RNA_Spearman_Rho_CytoscapeInput.txt");

$line = <IN>; # Get Time Point Name Column.
chomp $line;
@colnames = split("\t", $line);
shift(@colnames);
print OUT "AA_Sample_1\tAA_Sample_2\tSpearman_Rho\n";

$i = 0;
while ($line = <IN>) {
	chomp $line;
	@linearray = split("\t", $line);
	my $rowname = shift(@linearray);
	for (my $j = $i + 1; $j < scalar @linearray; ++$j) {
		if ($linearray[$j] >= 0.4 || $linearray[$j] <= -0.4) {
			print OUT "$rowname\t$colnames[$j]\t$linearray[$j]\n";
		}
	}
	++$i;
}

close IN;
close OUT;

exit;
