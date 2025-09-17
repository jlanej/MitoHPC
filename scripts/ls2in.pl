#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	# define variables
	my %opt;
	$opt{out}="out";

	my $result = GetOptions(
		"out=s"	=>	\$opt{out}
	);
        die "ERROR: $! " if (!$result);

	while(<>)
	{
		chomp;
		next unless(/\.bam$/ or /\.cram$/);
		my @F=split;

		# Extract sample name from BAM/CRAM header using samtools
		my $filepath = $F[-1];
		my $sample_name = `samtools samples "$filepath" 2>/dev/null | cut -f1 | head -n1`;
		chomp($sample_name);
		
		# Fallback to filename-based extraction if samtools fails
		if (!$sample_name) {
			if ($filepath =~ /.+\/(\S+)\./ or $filepath =~ /(\S+)\./) {
				$sample_name = $1;
			}
		}
		
		if ($sample_name) {
			print "$sample_name\t$filepath\t$opt{out}/$sample_name/$sample_name\n";
		}
	}
	exit 0;
}

