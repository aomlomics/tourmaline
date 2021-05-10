#!/usr/bin/perl

my $usage = qq{

    fastaLengths.pl <infile>

        <infile>: multi-FASTA file with or without line breaks

        This script takes a multi-fasta file with or without line breaks and
        prints the sequence header (up to the first space) and number of bases
        in each sequence, tab-delimited (.tsv), to standard output.


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

    # Store header line. If previous length exists, print previous length and
    # next header, else print (first) header. Initialize length.
    } elsif ($line =~ /^>[a-zA-Z0-9]*/) {
        $header = $line =~ s/>([a-zA-Z0-9]*).*\n/$1/r;
        if ($length) {
            print "$length\n";
            print "$header\t";
        } else {
            print "$header\t";
        }
        $length = 0;
        next;

    # Keep sequence line, add bases to length
    } else {
        chomp $line;
        $length += length($line);
    }
}

# Print last fasta length
print "$length\n";

exit;
