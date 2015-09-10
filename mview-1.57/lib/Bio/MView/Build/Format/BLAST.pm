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

# $Id: BLAST.pm,v 1.11 2005/12/12 20:42:48 brown Exp $

###########################################################################
#
# generic BLAST material
#
###########################################################################


###########################################################################
package Bio::MView::Build::Format::BLAST;

use vars qw(@ISA %Known_Parameter);
use NPB::Parse::Regexps;
use Bio::MView::Build::Search;
use strict;

@ISA = qw(Bio::MView::Build::Search);

#return the name of the underlying NPB::Parse::Stream parser
sub parser { 'BLAST' }

%Known_Parameter = 
    (
     #name        => [ format,     	   default  ]
						   
     #BLAST* display various HSP selections
     'hsp'        => [ '\S+',              'ranked' ],

     #BLAST* (version 1)					   
     'maxpval'    => [ $RX_Ureal,  	   undef    ],
     'minscore'   => [ '\d+',      	   undef    ],
						   
     #BLAST* (version 2)					   
     'maxeval'    => [ $RX_Ureal,  	   undef    ],
     'minbits'    => [ '\d+',      	   undef    ],
     'cycle'      => [ [],         	   undef    ],
	
     #BLASTN (version 1, version 2)
     #BLASTX (version 2)
     'strand'     => [ [],         	   undef    ],
	
     #BLASTX/TBLASTX (version 1)
##   'frame'      => [ [],         	   undef    ],
    );

sub initialise_parameters {
    my $self = shift;
    $self->SUPER::initialise_parameters;
    $self->SUPER::initialise_parameters(\%Known_Parameter);
    $self->reset_cycle;
    $self->reset_strand;
##  $self->reset_frame;
}

sub set_parameters {
    my $self = shift;
    $self->SUPER::set_parameters(@_);
    $self->SUPER::set_parameters(\%Known_Parameter, @_);
    $self->reset_cycle;
    $self->reset_strand;
##  $self->reset_frame;
}

sub subheader {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s    if $quiet;
    if ($self->{'hsp'} eq 'all') {
	$s .= "HSP processing: all\n";
    } elsif ($self->{'hsp'} eq 'discrete') {
	$s .= "HSP processing: discrete\n";
    } else {
	$s .= "HSP processing: ranked\n";
    }
    $s;    
}

sub new {
    shift;    #discard type
    my $self = new Bio::MView::Build::Search(@_);
    my ($type, $p, $v, $file);

    #determine the real type from the underlying parser
    ($p, $v) = (lc $self->{'entry'}->{'format'},$self->{'entry'}->{'version'});

    $type = "Bio::MView::Build::Format::BLAST$v";
    ($file = $type) =~ s/::/\//g;
    require "$file.pm";
    
    $type .= "::$p";
    bless $self, $type;

    $self->initialise;
}

#initialise parse iteration scheduler variable(s). just do them all at once
#and don't bother overriding with specific methods. likewise the scheduler
#routines can all be defined here.
sub initialise {
    my $self = shift;

    #NCBI BLAST2/PSI-BLAST search cycle
    $self->{'do_cycle'}    = undef;    #required list of cycle numbers
    $self->{'cycle_idx'}   = undef;    #current index into 'do_cycle'
    $self->{'cycle_ptr'}   = undef;    #current cycle parse object ref

    #BLASTN strand orientations
    $self->{'strand_list'} = [ qw(+ -) ];    #strand orientations
    $self->{'do_strand'}   = undef;    #list of required strand
    $self->{'strand_idx'}  = undef;    #current index into 'do_strand'

##  #BLASTX reading frame
##  $self->{'frame_list'}  = [ qw(+1 +2 +3 -1 -2 -3) ];    #reading frames
##  $self->{'do_frame'}    = undef;    #list of required frames
##  $self->{'frame_idx'}   = undef;    #current index into 'do_frame'

    $self->initialise_parameters;      #other parameters done last

    $self;
}

sub cycle    { $_[0]->{'do_cycle'}->[$_[0]->{'cycle_idx'}-1] }
sub strand   { $_[0]->{'do_strand'}->[$_[0]->{'strand_idx'}-1] }
##sub frame    { $_[0]->{'do_frame'}->[$_[0]->{'frame_idx'}-1] }

sub reset_cycle {
    my $self = shift;

    #initialise scheduler loops and loop counters
    if (@{$self->{'cycle'}} < 1) {
	#empty list  - do all cycles
	$self->{'do_cycle'} = [ 1..$self->{'entry'}->count(qw(SEARCH)) ];
    } elsif ($self->{'cycle'}->[0] != 0) {
	#explicit cycle range
	$self->{'do_cycle'} = [ @{$self->{'cycle'}} ];
    } else {
	#default to last cycle
	$self->{'do_cycle'} = [ $self->{'entry'}->count(qw(SEARCH)) ];
    }

    if (defined $self->{'cycle_ptr'}) {
	#flag previous cycle parse for garbage collection
	$self->{'cycle_ptr'}->free;
	$self->{'cycle_ptr'} = undef;
    }
}

sub reset_strand {
    my $self = shift;

    #initialise scheduler loops and loop counters
    if (@{$self->{'strand'}} < 1 or $self->{'strand'}->[0] eq '*') {
	#empty list  - do all strand
	$self->{'do_strand'} = [ @{$self->{'strand_list'}} ];
    } else {
	#explicit strand range
	$self->{'do_strand'} = [ @{$self->{'strand'}} ];
    }
}

##sub reset_frame {
##    my $self = shift;
##
##    #initialise scheduler loops and loop counters
##    if (@{$self->{'frame'}} < 1 or $self->{'frame'}->[0] eq '*') {
##	  #empty list  - do all frames
##	  $self->{'do_frame'} = [ @{$self->{'frame_list'}} ];
##    } else {
##	  #explicit frame range
##	  $self->{'do_frame'} = [ @{$self->{'frame'}} ];
##    }
##}

sub next_cycle {
    my $self = shift;

    #first pass?
    $self->{'cycle_idx'} = 0    unless defined $self->{'cycle_idx'};
    
    #normal pass: post-increment cycle counter
    if ($self->{'cycle_idx'} < @{$self->{'do_cycle'}}) {
	return $self->{'do_cycle'}->[$self->{'cycle_idx'}++];
    }

    #finished loop
    $self->{'cycle_idx'} = undef;
}

sub next_strand {
    my $self = shift;

    #first pass?
    $self->{'strand_idx'} = 0    unless defined $self->{'strand_idx'};
    
    #normal pass: post-increment strand counter
    if ($self->{'strand_idx'} < @{$self->{'do_strand'}}) {
	return $self->{'do_strand'}->[$self->{'strand_idx'}++];
    }

    #finished loop
    $self->{'strand_idx'} = undef;
}

##sub next_frame {
##    my $self = shift;
##
##    #first pass?
##    $self->{'frame_idx'} = 0    unless defined $self->{'frame_idx'};
##    
##    #normal pass: post-increment frame counter
##    if ($self->{'frame_idx'} < @{$self->{'do_frame'}}) {
##	return $self->{'do_frame'}->[$self->{'frame_idx'}++];
##    }
##
##    #finished loop
##    $self->{'frame_idx'} = undef;
##}

sub schedule_by_cycle {
    my ($self, $next) = shift;
    if (defined ($next = $self->next_cycle)) {
	return $next;
    }
    return undef;           #tell parser    
}

sub schedule_by_strand {
    my ($self, $next) = shift;
    if (defined ($next = $self->next_strand)) {
	return $next;
    }
    return undef;           #tell parser
}

##sub schedule_by_frame {
##    my ($self, $next) = shift;
##    if (defined ($next = $self->next_frame)) {
##	return $next;
##    }
##    return undef;           #tell parser
##}

sub schedule_by_cycle_and_strand {
    my ($self, $next) = (@_, 1);

    if (defined $self->{'cycle_idx'}) {
	#keep current cycle
	if (! defined $self->{'strand_idx'}) {
	    #strands finished: goto next cycle
	    $next = $self->next_cycle;
	}
    } else {
	#goto first cycle
	$next = $self->next_cycle;
    }
    
    #test the cycle
    if (defined $next) {
	#current/new cycle: goto next/first strand
	$next = $self->next_strand;
    } else {
	#all cycles finished
	return undef;           #tell parser    
    }

    #test the strand
    if (defined $next) {
	#ready to parse
	return $next;
    }

    #all cycles finished
    return undef;           	#tell parser    
}

##sub schedule_by_cycle_and_frame {
##    my ($self, $next) = (@_, 1);
##
##    if (defined $self->{'cycle_idx'}) {
##	#keep current cycle
##	if (! defined $self->{'frame_idx'}) {
##	    #frames finished: goto next cycle
##	    $next = $self->next_cycle;
##	}
##    } else {
##	#goto first cycle
##	$next = $self->next_cycle;
##    }
##    
##    #test the cycle
##    if (defined $next) {
##	#current/new cycle: goto next/first frame
##	$next = $self->next_frame;
##    } else {
##	#all cycles finished
##	return undef;           #tell parser    
##    }
##
##    #test the frame
##    if (defined $next) {
##	#ready to parse
##	return $next;
##    }
##
##    #all cycles finished
##    return undef;           	#tell parser    
##}

#override base class method to process query row differently
sub build_rows {
    my $self = shift;
    my ($lo, $hi, $i);

    #first, compute alignment length from query sequence in row[0]
    ($lo, $hi) = $self->set_range($self->{'index2row'}->[0]);
    
    #warn "range ($lo, $hi)\n";
       
    #query row contains missing query sequence, rather than gaps
    $self->{'index2row'}->[0]->assemble($lo, $hi, 'X');

    #assemble sparse sequence strings for all rows
    for ($i=1; $i < @{$self->{'index2row'}}; $i++) {
	$self->{'index2row'}->[$i]->assemble($lo, $hi, $self->{'gap'});
    }
    $self;
}


###########################################################################
package Bio::MView::Build::Row::BLAST;

use Bio::MView::Build::Row;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row);

sub posn1 {
    my $qfm = $_[0]->{'seq'}->fromlabel1;
    my $qto = $_[0]->{'seq'}->tolabel1;
    return "$qfm:$qto";
}

sub posn2 {
    my $hfm = $_[0]->{'seq'}->fromlabel2;
    my $hto = $_[0]->{'seq'}->tolabel2;
    return "$hfm:$hto"    if defined $_[0]->num and $_[0]->num;
    return '';
}

#fragment sort(worst to best): 1) increasing score, 2) increasing length
sub sort {
    $_[0]->{'frag'} =
	[
	 sort {
	     my $c = $a->[7] <=> $b->[7];                   #compare score
	     return $c    if $c != 0;
	     return length($a->[0]) <=> length($b->[0]);    #compare length
	 } @{$_[0]->{'frag'}}
	];
    $_[0];
}

sub assemble_blastp {
    my $self = shift;
    my ($i, $tmp);

    #query:     protein
    #database:  protein
    #alignment: protein x protein
    #query numbered in protein units
    #sbjct numbered in protein units
    #query orientation: +
    #sbjct orientation: +

    #processing steps:
    #  (1) assemble frags
   
    $self->SUPER::assemble(@_);
    $self;
}

sub assemble_blastn {
    my $self = shift;
    my ($i, $tmp);

    #query:     dna
    #database:  dna
    #alignment: dna x dna
    #query numbered in dna units
    #sbjct numbered in dna units
    #query orientation: +-
    #sbjct orientation: +-

    #processing steps:
    #if query -
    #  (1) reverse assembly position numbering
    #  (2) reverse each frag
    #  (3) assemble frags
    #  (4) reverse assembly
    #if query +
    #  (1) assemble frags
    
    if ($self->{'query_orient'} =~ /^\-/) {
	#stage (1,2,3,4)
	$self->SUPER::assemble(@_, 1);
    } else {
	#stage (1)
	$self->SUPER::assemble(@_);
    }
    $self;
}

sub assemble_blastx {
    my $self = shift;
    my ($i, $tmp);

    #query:     dna
    #database:  protein
    #alignment: protein x protein
    #query numbered in dna units
    #sbjct numbered in protein units
    #query orientation: +-
    #sbjct orientation: +

    #processing steps:
    #if query -
    #  (1) convert to protein units
    #  (2) reverse assembly position numbering
    #  (3) reverse each frag
    #  (4) assemble frags
    #  (5) reverse assembly
    #if query +
    #  (1) convert to protein units
    #  (2) assemble frags
    
    if ($self->{'query_orient'} =~ /^\-/) {
	#stage (1)
	for ($i=0; $i < @{$self->{'frag'}}; $i++) {
	    #start' = int((start+2)/3); stop' = int(stop/3)
	    $self->{'frag'}->[$i]->[2] = int(($self->{'frag'}->[$i]->[2]+2)/3);
	    $self->{'frag'}->[$i]->[1] = int($self->{'frag'}->[$i]->[1]/3);
	}
	#stage (2,3,4,5)
	$self->SUPER::assemble(@_, 1);
    } else {
	#stage (1)
	for ($i=0; $i < @{$self->{'frag'}}; $i++) {
	    #start' = int((start+2)/3); stop' = int(stop/3)
	    $self->{'frag'}->[$i]->[1] = int(($self->{'frag'}->[$i]->[1]+2)/3);
	    $self->{'frag'}->[$i]->[2] = int($self->{'frag'}->[$i]->[2]/3);
	}
	#stage (2)
	$self->SUPER::assemble(@_);
    }
    $self;
}

sub assemble_tblastn {
    my $self = shift;
    my ($i, $tmp);

    #query:     protein
    #database:  dna
    #alignment: protein x protein
    #query numbered in protein units
    #sbjct numbered in dna units
    #query orientation: +
    #sbjct orientation: +-

    #processing steps:
    #  (1) assemble frags
    
    $self->SUPER::assemble(@_);
    $self;
}

sub assemble_tblastx {
    my $self = shift;
    my ($i, $tmp);

    #query:     dna
    #database:  dna
    #alignment: protein x protein
    #query numbered in dna units
    #sbjct numbered in dna units
    #query orientation: +-
    #sbjct orientation: +-

    #processing steps:
    #if query -
    #  (1) convert to protein units
    #  (2) reverse assembly position numbering
    #  (3) reverse each frag
    #  (4) assemble frags
    #  (5) reverse assembly
    #if query +
    #  (1) convert to protein units
    #  (2) assemble frags
    
    if ($self->{'query_orient'} =~ /^\-/) {
	#stage (1)
	for ($i=0; $i < @{$self->{'frag'}}; $i++) {
	    #start' = int((start+2)/3); stop' = int(stop/3)
	    $self->{'frag'}->[$i]->[2] = int(($self->{'frag'}->[$i]->[2]+2)/3);
	    $self->{'frag'}->[$i]->[1] = int($self->{'frag'}->[$i]->[1]/3);
	}
	#stage (2,3,4,5)
	$self->SUPER::assemble(@_, 1);
    } else {
	#stage (1)
	for ($i=0; $i < @{$self->{'frag'}}; $i++) {
	    #start' = int((start+2)/3); stop' = int(stop/3)
	    $self->{'frag'}->[$i]->[1] = int(($self->{'frag'}->[$i]->[1]+2)/3);
	    $self->{'frag'}->[$i]->[2] = int($self->{'frag'}->[$i]->[2]/3);
	}
	#stage (2)
	$self->SUPER::assemble(@_);
    }
    $self;
}


###########################################################################
1;
