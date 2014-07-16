#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Time::Piece;
use Time::Seconds;
use XML::Simple;
use List::Util qw(min max sum);
use Scalar::Util qw(openhandle);
use Bio::Root::Version;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use FindBin;

my $bls_name = "proteins_1.bls";
my $bls = Bio::SearchIO->new(-file=>$bls_name, -format=>"blast");
while (my $res = $bls->next_result) {
    my $hit = $res->next_hit or next;
    my($pid,$prod,$gene,$EC) = ($res->query_name, $hit->description, '', '');
    my $hsp_count = 0;
    while(my $hsp = $hit->next_hsp()){
	$hsp_count++;
	print $hsp_count, "\t", $pid, "\t", $hit->name, "\t", $hsp->percent_identity, "\t", $hsp->evalue , "\t", $prod, "\n"
    }
}
