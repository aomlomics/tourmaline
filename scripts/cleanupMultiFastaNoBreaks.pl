#!/usr/bin/perl

use strict;
use warnings;

my $usage = qq{

    cleanupMultiFastaNoBreaks.pl <infile>

        <infile>: multi-FASTA file

        This script takes a multi-fasta file, removes header after the
        first non-word character, removes line breaks from sequences, and
        prints to standard output.

	WARNING: This version will print headers lacking sequences.

        Last edit: August 2018


};

die $usage unless @ARGV >= 1;

my $filename;
my @fasta;
my $flag = 0;
my $header;
my $sequence;
my $newfasta;

# Open file, read data, close
($filename) = @ARGV;
chomp $filename;
unless (open (FILE, $filename)) {
    print "The file $filename does not seem to exist!\n";
}
@fasta = <FILE>;
close FILE;

foreach my $line (@fasta) {

    # Discard blank line
    if ($line =~ /^\s*$/) {
        next;

    # Discard comment line
    } elsif ($line =~ /^\s*#/) {
        next;

    # Keep line, DO add pre-newline, store header up to first non-word character, append to $newfasta
    } elsif ($line =~ /^>\w*/) {
        if ($flag == 1) {$newfasta .= "\n"}
        $header = "$&" . "\n";
        $newfasta .= $header;
        next;

    # Keep line, store as sequence string, DO remove whitespace, append to $newfasta
    } else {
        $sequence = $line;
        $sequence =~ s/\s//g;
        $newfasta .= $sequence;
        $flag = 1;
    }
}

# DO add terminal newline
$newfasta .= "\n";

print $newfasta;

exit;
