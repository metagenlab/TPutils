#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Std qw'getopts';
use Data::Dumper;
use Carp;

use 5.010;

use Bio::TreeIO;

my %opts;
getopts( "i:o:", \%opts ) or print_usage();


(my $new_name = $opts{i}) =~ s/\..*$/\.phyloxml/;

print "Converting $opts{i}\n";

my $ti = Bio::TreeIO->new(-format => 'newick', -file => $opts{i});

my $tree = $ti->next_tree;

print "Writing to $new_name\n";

my $to = Bio::TreeIO->new(-format => 'phyloxml', -file => ">$new_name");

$to->write_tree($tree);
