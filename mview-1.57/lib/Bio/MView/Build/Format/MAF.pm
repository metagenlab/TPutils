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

# $Id: MAF.pm,v 1.1 2013/12/01 18:56:16 npb Exp $

###########################################################################
package Bio::MView::Build::Row::MAF;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row);

sub new {
    my $type = shift;
    my ($num, $id, $desc, $seq, $start, $size, $strand, $srcsize) = @_;
    my $self = new Bio::MView::Build::Row($num, $id, $desc, $seq);
    $self->{'start'}   = $start;
    $self->{'size'}    = $size;
    $self->{'strand'}  = $strand;
    $self->{'srcsize'} = $srcsize;
    bless $self, $type;
}

sub data {
    return sprintf("%8s %8s %6s %10s",
		   $_[0]->{'start'}, $_[0]->{'size'}, $_[0]->{'strand'},
		   $_[0]->{'srcsize'})
	if $_[0]->num;
    return sprintf("%8s %8s %6s %10s", 'start', 'size', 'strand', 'srcsize');
}

sub rdb {
    my ($self, $mode) = (@_, 'data');
    my $s = $self->SUPER::rdb($mode);
    return join("\t", $s, $self->{'start'}, $self->{'size'}, $self->{'strand'},
		$self->{'srcsize'})
	if $mode eq 'data';
    return join ("\t", $s, 'start', 'size', 'strand', 'srcsize')
	if $mode eq 'attr';
    return join ("\t", $s, '8N', '8N', '6S', '10S')
	if $mode eq 'form';
    return '';
}


###########################################################################
package Bio::MView::Build::Format::MAF;

use vars qw(@ISA %Known_Parameter);
use Bio::MView::Build::Align;
use Bio::MView::Build::Row;
use strict;

@ISA = qw(Bio::MView::Build::Align);

#the name of the underlying NPB::Parse::Format parser
sub parser { 'MAF' }

my %Known_Parameter = 
    (
     #name        => [ format,               default ]
     'block'      => [ [],                   undef   ],
    );

sub initialise_parameters {
    my $self = shift;
    $self->SUPER::initialise_parameters;
    $self->SUPER::initialise_parameters(\%Known_Parameter);

    $self->reset_block;
}

sub set_parameters {
    my $self = shift;
    $self->SUPER::set_parameters(@_);
    $self->SUPER::set_parameters(\%Known_Parameter, @_);

    $self->reset_block;
}

sub new {
    my $type = shift;
    my $self = new Bio::MView::Build::Align(@_);

    #MAF ordinal block number: counted 1..N whereas the actual
    #block has its own 'number' field which is reported by subheader().
    $self->{'do_block'}  = undef;    #required list of block numbers
    $self->{'block_idx'} = undef;    #current index into 'do_block'
    $self->{'block_ptr'} = undef;    #current block parse object ref

    bless $self, $type;
}

sub block   { $_[0]->{'do_block'}->[$_[0]->{'block_idx'}-1] }

sub reset_block {
    my $self = shift;

    #initialise scheduler loops and loop counters
    if (! defined $self->{'do_block'}) {
	if (@{$self->{'block'}} < 1) {
	    #empty list  - do all blocks
	    $self->{'do_block'} = [ 1..$self->{'entry'}->count(qw(BLOCK)) ];
	} elsif ($self->{'block'}->[0] != 0) {
	    #explicit block range
	    $self->{'do_block'} = [ @{$self->{'block'}} ];
	} else {
	    #default to first block
	    $self->{'do_block'} = [ 1 ];
	}
    } else {
	#flag previous block parse for garbage collection
	$self->{'block_ptr'}->free;
	$self->{'block_ptr'} = undef;
    }
}

sub next_block {
    my $self = shift;

    #first pass?
    $self->{'block_idx'} = 0    unless defined $self->{'block_idx'};
    
    #normal pass: post-increment block counter
    if ($self->{'block_idx'} < @{$self->{'do_block'}}) {
	return $self->{'do_block'}->[$self->{'block_idx'}++];
    }

    #finished loop
    $self->{'block_idx'} = undef;
}

sub schedule_by_block {
    my ($self, $next) = shift;

    if (defined ($next = $self->next_block)) {
	return $next;
    }
    return undef;           #tell parser    
}

sub subheader {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s    if $quiet;
    $s .= "Block: " . $self->block . " (score:  $self->{'block_ptr'}->{'score'})\n";
    $s;    
}

sub parse {
    my $self = shift;
    my ($rank, $use, $desc, $seq, @hit) = (0);

    return  unless defined $self->schedule_by_block;

    $self->{'block_ptr'} = $self->{'entry'}->parse("BLOCK[@{[$self->block]}]");

    #block doesn't exist?
    return unless defined $self->{'block_ptr'};

    foreach my $row (@{$self->{'block_ptr'}->{'row'}}) {

	$rank++;

	#check row wanted, by rank OR identifier OR row count limit
	last  if ($use = $self->use_row($rank, $rank, $row->{'id'})) < 0;
	next  unless $use;

	#warn "KEEP: ($rank,$row->{'id'})\n";

	push @hit, new Bio::MView::Build::Row::MAF($rank,
						   $row->{'id'},
						   '',
						   $row->{'seq'},
						   $row->{'start'},
						   $row->{'size'},
						   $row->{'strand'},
						   $row->{'srcsize'},
	    );
    }

    #map { $_->print } @hit;
    
    return \@hit;
}


###########################################################################
1;
