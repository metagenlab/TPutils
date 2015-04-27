# -*- perl -*-
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

# $Id: ssearch.pm,v 1.8 2013/09/09 21:31:05 npb Exp $

###########################################################################
package NPB::Parse::Format::FASTA2::ssearch;

use NPB::Parse::Format::FASTA2;
use strict;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA2);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package NPB::Parse::Format::FASTA2::ssearch::HEADER;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA2::HEADER);

sub new {
    my $self = shift;
    $self = $self->SUPER::new(@_);
    #assume the query identifier is the same as the query filename
    if ($self->{query} eq '' and $self->{queryfile} ne '') {
	$self->{query} = $self->{queryfile};
    }
    return $self;
}

###########################################################################
package NPB::Parse::Format::FASTA2::ssearch::RANK;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA::RANK);

sub new {
    my $type = shift;
    if (@_ < 2) {
	#at least two args, ($offset, $bytes are optional).
	NPB::Message::die($type, "new() invalid arguments (@_)");
    }
    my ($parent, $text, $offset, $bytes) = (@_, -1, -1);
    my ($self, $line, $record);
    
    $self = new NPB::Parse::Record($type, $parent, $text, $offset, $bytes);
    $text = new NPB::Parse::Record_Stream($self);

    #ranked search hits
    while (defined ($line = $text->next_line)) {
	
	next    if $line =~ /$NPB::Parse::Format::FASTA2::RANK_START/o;

	if($line =~ /^
	   \s*
	   ([^\s]+)                #id
	   \s+
	   (.*)                    #description: possibly empty
	   \s+
	   \(\s*(\S+)\)            #hit sequence length
           \s*
	   (?:\[(\S)\])?           #frame
	   \s*
	   (\S+)                   #s-w score
	   \s+
	   (\S+)                   #Z-score
	   \s+
	   (\S+)                   #E-value
	   \s*
	   $/xo) {
	    
	    $self->test_args($line, $1, $3, $5,$6,$7);
	    
	    push(@{$self->{'hit'}},
		 { 
		  'id'      => NPB::Parse::Record::clean_identifier($1),
		  'desc'    => $2,
		  #ignore $3
		  'frame'   => NPB::Parse::Format::FASTA::parse_frame($4),
		  'orient'  => NPB::Parse::Format::FASTA::parse_orient($4),
		  'score'   => $5,
		  'zscore'  => $6,
		  'expect'  => $7,
		 });
	    next;
	}
    
	#blank line or empty record: ignore
	next    if $line =~ /$NPB::Parse::Format::FASTA2::NULL/o;
	
	#default
	$self->warn("unknown field: $line");
    }
    $self;
}


###########################################################################
package NPB::Parse::Format::FASTA2::ssearch::TRAILER;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA::TRAILER);


###########################################################################
package NPB::Parse::Format::FASTA2::ssearch::MATCH;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA::MATCH);


###########################################################################
package NPB::Parse::Format::FASTA2::ssearch::MATCH::SUM;

use vars qw(@ISA);
use NPB::Parse::Regexps;

@ISA   = qw(NPB::Parse::Format::FASTA::MATCH::SUM);

sub new {
    my $type = shift;
    if (@_ < 2) {
	#at least two args, ($offset, $bytes are optional).
	NPB::Message::die($type, "new() invalid arguments (@_)");
    }
    my ($parent, $text, $offset, $bytes) = (@_, -1, -1);
    my ($self, $line, $record);
    
    $self = new NPB::Parse::Record($type, $parent, $text, $offset, $bytes);
    $text = new NPB::Parse::Record_Stream($self);

    $line = $text->next_line;
	
    if ($line =~ /^
	>*
	(\S+)                      #id
	\s+
	(.*)                       #description: possibly empty
	\s+
	\(\s*(\d+)\s*(?:aa|nt)\)   #length
	\s*
	$/xo) {

	$self->test_args($line, $1, $3);

	(
	 $self->{'id'},
	 $self->{'desc'},
	 $self->{'length'},
	) = (NPB::Parse::Record::clean_identifier($1),
	     NPB::Parse::Record::strip_english_newlines($2), $3);
    } else {
	$self->warn("unknown field: $line");
    }

    $line = $text->next_line;
    
    if ($line =~ /^
	\s*
	z-score\:\s*(\S+)      #z
	\s*
	Expect\:\s*(\S+)       #E
	\s*
	$/xo) {

	$self->test_args($line,$1,$2);

	(
	 $self->{'zscore'},
	 $self->{'expect'},
	) = ($1,$2);
    } else {
	$self->warn("unknown field: $line");
    }
    
    $line = $text->next_line;

    if ($line =~ /^
	(?:Smith-Waterman\s+score:\s*(\d+);)?    #sw score
	\s*($RX_Ureal)%                          #percent identity
	\s*identity\s+in\s+(\d+)                 #overlap length
	\s+(?:aa|nt)\s+overlap
	\s*
	$/xo) {

	$self->test_args($line,$2,$3);
	
	(
	 $self->{'score'},
	 $self->{'id_percent'},
	 $self->{'overlap'},
	) = ((defined $1 ? $1 : 0), $2, $3);
    } else {
	$self->warn("unknown field: $line");
    }
    
    $self;
}


###########################################################################
package NPB::Parse::Format::FASTA2::ssearch::MATCH::ALN;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA2::MATCH::ALN);


###########################################################################
1;
