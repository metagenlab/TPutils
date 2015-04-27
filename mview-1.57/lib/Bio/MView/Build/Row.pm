###########################################################################
#
# Copyright (C) 1997-2015 Nigel P. Brown
# 
# (i) License
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
# 
# (ii) Contacts
# 
#  Project Admin:      Nigel P. Brown
#  Email:              biomview@gmail.com
#  Project URL:        http://bio-mview.sourceforge.net
# 
# (iii) Citation
# 
#  Please acknowledge use of this Program by citing the following reference in
#  any published work including web-sites:
#  
#   Brown, N.P., Leroy C., Sander C. (1998) MView: A Web compatible database
#   search or multiple alignment viewer. Bioinformatics. 14(4):380-381.
#  
#  and provide a link to the MView project URL given above under 'Contacts'.
#
###########################################################################

# $Id: Row.pm,v 1.14 2013/09/09 21:31:04 npb Exp $

###########################################################################
package Bio::MView::Build::Row;

use vars qw($Default_IdWidth $Default_TextWidth);
use Bio::MView::Sequence;
use strict;

$Default_IdWidth   = 30;    #default width to truncate 'id' field
$Default_TextWidth = 30;    #default width to truncate 'text' field

sub new {
    my $type = shift;
    my ($num, $id, $desc, $seq) = (@_, undef);
    my $self = {};

    bless $self, $type;

    #strip non-identifier leading rubbish:  >  or /:
    $id =~ s/^(>|\/:)//;

    $self->{'rid'}  = $id;                      #supplied identifier
    $self->{'uid'}  = $self->uniqid($num, $id); #unique compound identifier
    
    #ensure identifier is non-null (for Build::map_id())
    $id = ' '  unless $id =~ /./;

    #set row 'subtype' information
    if ($id =~ /^\#/) {
	$self->{'type'} = 'special';  #leading hash: user supplied row
	$id =~ s/^\#//;               #and strip it
    } else {
	$self->{'type'} = undef;                    
    }

    $self->{'num'}  = $num;                     #row number/string
    $self->{'cid'}  = $id;                      #cleaned identifier
    $self->{'desc'} = $desc;                    #description string
    $self->{'frag'} = [];                       #list of fragments

    $self->{'seq'}  = new Bio::MView::Sequence; #finished sequence

    $self->add_frag($seq)    if defined $seq;

    $self->{'url'}  = Bio::SRS::srsLink($self->{'cid'});  #url

    $self;
}

#sub DESTROY { warn "DESTROY $_[0]\n" }

#methods returning standard strings for use in generic output modes
sub rid  { $_[0]->{'rid'} }
sub uid  { $_[0]->{'uid'} }
sub cid  { $_[0]->{'cid'} }
sub num  { $_[0]->{'num'} }
sub url  { $_[0]->{'url'} }
sub sob  { $_[0]->{'seq'} }

sub seq  {
    return $_[0]->{'seq'}->string    if defined $_[0]->{'seq'};
    return '';
}

sub desc { $_[0]->{'desc'} }

sub data { '' }

sub text {
    my $w = defined $_[1] ? $_[1] : $Default_TextWidth;
    $w = length $_[0]->{'desc'}    if $w > length $_[0]->{'desc'};
    sprintf("%-${w}s", $_[0]->truncate($_[0]->{'desc'}, $w));
}

sub posn1 { '' }
sub posn2 { '' }

sub uniqid { "$_[1]\034/$_[2]" }

sub print {
    sub _format {
	my ($self, $k, $v) = @_;
	$v = 'undef' unless defined $v;
	$v = "'$v'" if $v =~ /^\s*$/;
	return sprintf("  %-15s => %s\n", $k, $v)
    }
    my $self = shift;
    warn "$self\n";
    map { warn $self->_format($_, $self->{$_}) } sort keys %{$self};
    $self;
}

sub truncate {
    my ($self, $s, $n, $t) = (@_, $Default_TextWidth);
    $t = substr($s, 0, $n);
    substr($t, -3, 3) = '...'    if length $s > $n;
    $t;
}

#routine to sort 'frag' list: default is null
sub sort {$_[0]}

#modify the extra 'type' information
sub set_subtype { $_[0]->{'type'} = $_[1] }

#add a sequence fragment to the 'frag' list with value and positions given
#by first three args. use default positions if called with one arg. other
#optional arguments are special to any subclass of Row.
sub add_frag {
    my $self = shift;
    my ($frag, $qry_from, $qry_to) = (shift, shift, shift);

    $qry_from = 1               unless defined $qry_from;
    $qry_to   = length $frag    unless defined $qry_to;

    push @{$self->{'frag'}}, [ \$frag, $qry_from, $qry_to, @_ ];

    #warn "@{$self->{'frag'}->[ $#{$self->{'frag'}} ]}\n";

    $self;
}

sub count_frag { scalar @{$_[0]->{'frag'}} }

#compute the maximal positional range of a row
sub range {
    my ($self, $orient) = (@_, '+');
    return $self->range_plus    if $orient =~ /^\+/;
    return $self->range_minus;
}

#compute the maximal positional range of a forward row (numbered forwards)
sub range_plus {
    my $self = shift;
    my ($lo, $hi, $frag);

    return (0, 0)    unless @{$self->{'frag'}};

    $lo = $hi = $self->{'frag'}->[0]->[1];

    foreach $frag (@{$self->{'frag'}}) {
	$lo = ($frag->[1] < $lo ? $frag->[1] : $lo);
	$hi = ($frag->[2] > $hi ? $frag->[2] : $hi);
    }
    #warn "range_plus ($lo, $hi)\n";
    ($lo, $hi);
}

#compute the maximal positional range of a reversed row (numbered backwards)
sub range_minus {
    my $self = shift;
    my ($lo, $hi, $frag);

    return (0, 0)    unless @{$self->{'frag'}};

    $lo = $hi = $self->{'frag'}->[0]->[1];

    foreach $frag (@{$self->{'frag'}}) {
	$lo = ($frag->[2] < $lo ? $frag->[2] : $lo);
	$hi = ($frag->[1] > $hi ? $frag->[1] : $hi);
    }
    #warn "range_minus ($lo, $hi)\n";
    ($lo, $hi);
}

sub assemble {
    my ($self, $lo, $hi, $gap, $reverse) = (@_, 0);
    $self->sort;                                   #fragment order
    $self->{'seq'}->reverse    if $reverse;        #before calling append() !
    $self->{'seq'}->append(@{$self->{'frag'}});    #assemble fragments
    $self->{'seq'}->set_range($lo, $hi);           #set sequence range
    $self->{'seq'}->set_pad($gap);
    $self->{'seq'}->set_gap($gap);
    $self;
}

sub set_pad { $_[0]->{'seq'}->set_pad($_[1]) }
sub set_gap { $_[0]->{'seq'}->set_gap($_[1]) }
sub set_spc { $_[0]->{'seq'}->set_spc($_[1]) }

sub rdb {
    my ($self, $mode) = (@_, 'data');
    if ($mode eq 'data') {
	return join "\t", $self->num, $self->cid, $self->url, $self->seq, $self->desc;
    }
    if ($mode eq 'attr') {
	return join "\t", 'num', 'cid', 'url', 'seq', 'desc';
    }
    if ($mode eq 'form') {
	return join "\t", '4N', '30S', '100S', '500S', '500S';
    }
    '';
}

sub pearson {
    my $self = shift;
    my ($s, $p, $d, $i) = ($self->seq);
    $p = ">";
    $d = $self->num;
    $p .= ((defined $d and $d ne '') ? "$d;" : "query;");
    $p .= $self->cid;
    $d = $self->desc;  $p .= " $d"  if $d ne '';
    $d = $self->data;  $p .= " $d"  if $d ne '';
    $d = $self->posn1; $p .= " $d"  if $d ne '';
    $d = $self->posn2; $p .= " $d"  if $d ne '';
    $p .= "\n";
    for ($i=0; $i<length($s); $i+=70) {
        $p .= substr($s, $i, 70) . "\n";
    }
    $p;
}

sub pir {
    my $self = shift;
    my ($s, $p, $q, $d, $i) = ($self->seq);
    $p = ">P1;";
    $d = $self->num;
    $p .= ((defined $d and $d ne '') ? "$d;" : "query;");
    $p .= $self->cid . "\n";
    $d = $self->desc;  $p .= " $d"  if $d ne '';
    $d = $self->data;  $p .= " $d"  if $d ne '';
    $d = $self->posn1; $p .= " $d"  if $d ne '';
    $d = $self->posn2; $p .= " $d"  if $d ne '';
    for ($i=0; $i<length($s); $i+=60) {
	$q .= "\n" . substr($s, $i, 60);
    }
    $q .= "\n"    if length($q) % 61 < 1 and $q ne '';
    $q .= "*\n";
    $p .= $q;
}

sub plain {
    my ($self, $w) = (@_, $Default_IdWidth);
    my ($s, $p) = ($self->seq, $self->cid);
    $p = substr($p, 0, $w)  if length $p > $w;
    $p = sprintf("%-${w}s", $p) . " $s\n";
}


###########################################################################
1;
