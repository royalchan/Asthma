#!/usr/bin/perl

use strict;
use warnings;

####################################
# This script is to combine the histone marker peaks.
####################################

use FileHandle;

my $line;
my @linearray;

open(FILES, "04_AA_RNA_Hits_Bed_HistoneMarkPeak_Combiner_Files.txt") || die;
my @filenames;
my @outnames;
while ($line = <FILES>) {
	chomp $line;
	push (@filenames, $line);
	my $tmp = $line;
	$tmp =~ s/\.bed//;
	$tmp .= "_peakcombined.bed";
	push (@outnames, $tmp);
}
close FILES;

my $fh = FileHandle -> new;
my $oh = FileHandle -> new;

for (my $i = 0; $i < (scalar @filenames); ++$i) {

	$fh -> open ("< $filenames[$i]");
	$oh -> open ("> $outnames[$i]");
	
	# Output header.
	$line = <$fh>;
	print $oh $line;
	
	# Final printed combined values.
	# Initiate these values with values from the first line.
	$line = <$fh>;
	chomp $line;
	@linearray = split("\t", $line);	
	my $chrprt = $linearray[0];
	my $startprt = $linearray[1];
	my $stopprt = $linearray[2];
	my $sumpk = $linearray[3];
	my $count = 1;
	
	while ($line = <$fh>) {
	
		chomp $line;	
		@linearray = split("\t", $line);
		
		if (($linearray[0] eq $chrprt) && ($linearray[1] eq $stopprt)) {
			$stopprt = $linearray[2];
			$sumpk += $linearray[3];
			++$count;
		} else {
			my $avepk = $sumpk / $count;
			print $oh "$chrprt\t$startprt\t$stopprt\t$avepk\n";
			$chrprt = $linearray[0];
			$startprt = $linearray[1];
			$stopprt = $linearray[2];
			$sumpk = $linearray[3];
			$count = 1;
		}
		
	}
	
	# Print last record.
	my $avepk = $sumpk / $count;
	print $oh "$chrprt\t$startprt\t$stopprt\t$avepk\n";

	$fh -> close;
	$oh -> close;
	
}

exit;
