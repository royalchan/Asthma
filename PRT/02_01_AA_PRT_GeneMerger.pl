#!/usr/bin/perl

use warnings;
use strict;

####################################
# This script is to merge the AA Proteome MS intensity by gene name.
# Created by Rui Chen, 12/10/2015.
####################################
my $line;
my @linearray;

open(IN, "AA_Allreplicate1_reloaded_12092015.txt") || die;
open(OUT, "> AA_Allreplicate1_reloaded_12092015_genecbd.txt") || die;
my $header = <IN>;
$header =~ s/^Accession\tProtein\ Name\t//;
print OUT $header;
$header = <IN>;
$header =~ s/^Accession\tProtein\ Name\t//;
print OUT $header;

my %geneint; # Hash of which keys are gene names and values are scalar of readings.
my %genecount; # Count for gene replicates of which keys are gene names and values are scalar of counts for each sample.

while ($line = <IN>) {
	
	chomp $line;
	@linearray = split("\t", $line);
	
	for (my $i = 3; $i < scalar @linearray; ++$i) {
		if (defined $geneint{$linearray[2]} || 0) {
			if ($linearray[$i] eq ".") {
				${$geneint{$linearray[2]}}[$i - 3] += 0;
				${$genecount{$linearray[2]}}[$i - 3] += 0;
			}	else {
				${$geneint{$linearray[2]}}[$i - 3] += $linearray[$i];
				${$genecount{$linearray[2]}}[$i - 3] += 1;
			}
		} else {
			if ($linearray[$i] eq ".") {
				${$geneint{$linearray[2]}}[$i - 3] = 0;
				${$genecount{$linearray[2]}}[$i - 3] = 0;
			} else {
				${$geneint{$linearray[2]}}[$i - 3] = $linearray[$i];
				${$genecount{$linearray[2]}}[$i - 3] = 1;
			}
		}
	}

}

# Now output file.
my @keys = sort {$a cmp $b} (keys %geneint);
foreach (@keys) {
	my @meanint; # Mean intensity;
	for (my $j = 0; $j < scalar @{$geneint{$_}}; ++$j) {
		if (${$genecount{$_}}[$j] != 0) {
			$meanint[$j] = ${$geneint{$_}}[$j] / ${$genecount{$_}}[$j];
		} else {
			$meanint[$j] = 0;
		}
	}
	my $reline = join("\t", @meanint);
	$reline = $_."\t".$reline."\n";
	print OUT $reline;		
}

close IN;
close OUT;

exit;
