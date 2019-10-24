#!/usr/bin/perl
# contiguous_fasta.pl -- splits fasta-formatted files into contiguous
# sequences of non-ambiguous bases
# Author: David Eccles (gringer) 2011 <david.eccles@mpi-muenster.mpg.de>

use warnings;
use strict;

my $id = "";
my $first = 1; # true
my @sequences = ();
my $seqNum = 0;

while(<>){
    chomp;
    if(substr($_,0,1) eq ">"){
        $id = $_;
        $seqNum = 0;
        print((($first)?"":"\n").$_."/".$seqNum++."\n");
        $first = 0; # false
    } else {
        @sequences = split(/N+/, $_);
        print(shift(@sequences)) unless !@sequences; # whole line could be N
        foreach my $sequence (@sequences){
            print("\n$id/".$seqNum++."\n$sequence");
        }
    }
}
print("\n");
