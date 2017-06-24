#!/usr/bin/env perl

## prepand copying information

$#ARGV >= 1 || die "usage: prepend-notice.pl <notice> <file 1> [<file 2> ...]\n";

$copying_src = shift @ARGV;

print "## applying notice in $copying_src\n";

$copying = `cat $copying_src`;

foreach (@ARGV) {
    print "## processing $_\n";
    $orig = `cat $_`;

    ## get rid of all CRs
    $orig =~ s/\r//gs;

    open(FOUT, ">$_");
    print FOUT "%% $_\n\n";
    print FOUT "$copying\n";
    print FOUT "$orig";
    close(FOUT);
}
