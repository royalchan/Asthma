#!/usr/bin/perl

use strict;
use warnings;

####################################
# This script is to merge the 5 AA proteome files.
# Created by Rui Chen, 11/01/2015.
# After this step use the unix command grep -v "\t\.\t" to get rid of the lines with missing values.
####################################

use FileHandle;

my $line;
my @linearray;

my @files = ("100515_PBMCAA1R1_Simplified.txt", "102615_PBMCAA2R1_Simplified.txt", "102215_PBMCAA6R1_Simplified.txt", "101015_PBMCAA8R1_Simplified.txt", "101715_PBMCAA8R1_Simplified.txt");
my $outfile = "AA_PRT_Combined_Results_20151101.txt";

my $fh = FileHandle -> new;
my $oh = FileHandle -> new;
$oh -> open("> $outfile");
my $header = "Accession";

# Define a hash who's keys are Accession Numbers and values being concatenated ratios.
my %cbd;

# Now process the files one by one.
for (my $i = 0; $i < scalar @files; ++$i) {

	$fh -> open("< $files[$i]");
	
	# Get and concatenate header.
	$line = <$fh>;
	chomp $line;
	$line =~ s/^Accession\t//;
	$header .= "\t$line";
	
	# Now process the data line by line.
	while ($line = <$fh>) {
		chomp $line;
		# Skip the empty lines;
		next if ($line =~ /\t\t/);
		# Get the key and the values, and populate the hash.
		@linearray = split("\t", $line);
		my $key = $linearray[0];
		shift @linearray;
		my $reline = join("\t", @linearray);
		if ($i == 0) {
			$cbd{$key} = "$reline";
		} elsif (defined $cbd{$key} || 0) {
			my @linearraycbd = split("\t", $cbd{$key});
			if (scalar @linearraycbd == 9 * $i) {
				$cbd{$key} .= "\t$reline";			
			} else {
				do {push(@linearraycbd, ".");} until (scalar @linearraycbd == 9 * $i);
				my $relinecbd = join("\t", @linearraycbd);
				$cbd{$key} = "$relinecbd\t$reline";
			}
		} else {
			my @linearraycbd = (".");	
			do {push(@linearraycbd, ".");} until (scalar @linearraycbd == 9 * $i);
			my $relinecbd = join("\t", @linearraycbd);
			$cbd{$key} = "$relinecbd\t$reline";
		}
	}
	
	$fh -> close;

}

# Now output the combined file.
# Print out header.
print $oh "$header\n";
# Print out sorted ratios by Accession Number.
my @keys = sort {$a cmp $b} (keys %cbd);
foreach (@keys) {
	# First fill in the accession number that misses values from the last file with dots.
	@linearray = split("\t", $cbd{$_});
	if (scalar @linearray < 45) {
		do {push(@linearray, ".");} until (scalar @linearray == 45);
	}
	my $reline = join("\t", @linearray);
	# Then print out the final line.
	print $oh "$_\t$reline\n";
}

$oh -> close;

exit;
