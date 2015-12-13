#!/usr/bin/perl

use warnings;
use strict;

####################################
# This script is to find Hannes Hits by gene name.
# Created by Rui Chen, 12/10/2015.
####################################
my $line;
my @linearray;

open(IN, "AA_Allreplicate1_reloaded_12092015_genecbd.txt") || die;
open(IN1, "genes_output_p_adj_anova_0_01.txt") || die;
open(IN2, "anova_alltwins_genes_output_p_simple_0_01.txt") || die;
open(OUT1, "> AA_Allreplicate1_reloaded_12092015_genecbd_p0_01_daonly.txt") || die;
open(OUT2, "> AA_Allreplicate1_reloaded_12092015_genecbd_p0_01_full.txt") || die;
my $header = <IN>;
print OUT1 $header;
print OUT2 $header;
$header = <IN>;
print OUT1 $header;
print OUT2 $header;

my %data;
while ($line = <IN>) {
	@linearray = split("\t", $line);
	$data{$linearray[0]} = $line;
}

# Now output file.
while ($line = <IN1>) {
	chomp $line;
	if (defined $data{$line} || 0) {print OUT1 $data{$line};}
}
while ($line = <IN2>) {
	chomp $line;
	if (defined $data{$line} || 0) {print OUT2 $data{$line};}
}

close IN;
close IN1;
close IN2;
close OUT1;
close OUT2;

exit;
