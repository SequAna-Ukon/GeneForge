#!/usr/bin/env perl
use strict;
use warnings;

# Check for input and output file arguments
if (@ARGV != 2) {
    die "Usage: $0 <input_file> <output_file>\n";
}

my ($input_file, $output_file) = @ARGV;

# Open input and output files
open my $in_fh, '<', $input_file or die "Cannot open input file $input_file: $!\n";
open my $out_fh, '>', $output_file or die "Cannot open output file $output_file: $!\n";

# Process file
my $line_count = 0;
while (my $line = <$in_fh>) {
    $line_count++;
    chomp $line;

    # Print first three lines (headers)
    if ($line_count <= 3) {
        print $out_fh "$line\n";
        next;
    }

    # Split line on tabs
    my @fields = split /\t/, $line;

    # Check if last field contains 'high' (case-insensitive)
    if (@fields && $fields[-1] =~ /high/i) {
        print $out_fh "$line\n";
    }
}

# Close files
close $in_fh;
close $out_fh;
