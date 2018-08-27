#!/usr/bin/perl

my $usage = qq{

    fastaLengthDist.pl <infile>

        <infile>: multi-FASTA file with or without line breaks

        This script takes a multi-fasta file with or without line breaks and
        prints the number of bases in each sequence to standard output.


};

die $usage unless @ARGV == 1;

$filename = $ARGV[0];

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

    # Discard header line but print previous length and initialize length
    } elsif ($line =~ /^>[a-zA-Z0-9]*/) {
	if ($length) { print "$length\n"; }
        $length = 0;
	next;

    # Keep sequence line, count bases, print
    } else {
        chomp $line;
	$length += length($line);
    }
}

#print last fasta length
print "$length\n";

exit;
