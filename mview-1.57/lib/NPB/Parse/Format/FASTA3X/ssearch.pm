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

# $Id: ssearch.pm,v 1.2 2013/09/09 21:31:06 npb Exp $

###########################################################################
package NPB::Parse::Format::FASTA3X::ssearch;

use NPB::Parse::Format::FASTA3X;
use strict;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3X);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package NPB::Parse::Format::FASTA3X::ssearch::HEADER;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA3X::HEADER);

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
package NPB::Parse::Format::FASTA3X::ssearch::RANK;

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
	
	next    if $line =~ /$NPB::Parse::Format::FASTA3X::RANK_START/o;

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
	   (\S+)                   #bits
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
		  'bits'    => $6,
		  'expect'  => $7,
		 });
	    next;
	}
    
	#blank line or empty record: ignore
	next    if $line =~ /$NPB::Parse::Format::FASTA3X::NULL/o;
	
	#default
	$self->warn("unknown field: $line");
    }
    $self;
}


###########################################################################
package NPB::Parse::Format::FASTA3X::ssearch::TRAILER;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA::TRAILER);


###########################################################################
package NPB::Parse::Format::FASTA3X::ssearch::MATCH;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA::MATCH);


###########################################################################
package NPB::Parse::Format::FASTA3X::ssearch::MATCH::SUM;

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

    my $lines = $text->scan_until_inclusive('^\s*Smith-Waterman');

    if ($lines =~ /^
	>*
	(\S+)                      	  #id
	\s+
	(.*)                       	  #description: possibly empty
	\s+
	\(\s*(\d+)\s*(?:aa|nt)\)   	  #length
	\s*
        (rev-comp)?                       #frame
        \s+
        s-w\s+opt:\s*(\d+)         	  #s-w opt
        \s+
        Z-score:\s*(\S+)           	  #Z-score
        \s+
        bits:\s*(\S+)              	  #bits
        \s+
        E\((?:\d+)?\):\s*(\S+)         	  #expect
        \s+
        Smith-Waterman\s+score:\s*(\S+);  #s-w score (same as opt?)
        \s+
        (\S+)%\s+identity                 #percent identity
        \s+
        \((\S+)%\s+(?:similar|ungapped)\) #percent similarity
        \s+in\s+
        (\d+)\s*(?:aa|nt)\s+overlap       #alignment overlap
        \s+
        \((\S+)\)                         #alignment ranges
        \s*
	$/xso) {

	$self->test_args($line, $1, $3, $5,$6,$7,$8,$9,$10,$11,$12,$13);

	(
	 $self->{'id'},
	 $self->{'desc'},
	 $self->{'length'},
	 $self->{'frame'},
	 $self->{'orient'},
	 $self->{'opt'},
	 $self->{'zscore'},
	 $self->{'bits'},
	 $self->{'expect'},
	 $self->{'score'},
	 $self->{'id_percent'},
	 $self->{'sim_percent'},
	 $self->{'overlap'},
	 $self->{'ranges'},
	) = (
	    NPB::Parse::Record::clean_identifier($1),
	    NPB::Parse::Record::strip_english_newlines($2),
	    $3,
	    NPB::Parse::Format::FASTA::parse_frame($4),
	    NPB::Parse::Format::FASTA::parse_orient($4),
	    $5, $6, $7, $8, $9, $10, $11, $12, $13,
	);
    } else {
	$self->warn("unknown field: $lines");
    }

    $self;
}


###########################################################################
package NPB::Parse::Format::FASTA3X::ssearch::MATCH::ALN;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA3X::MATCH::ALN);


###########################################################################
1;
