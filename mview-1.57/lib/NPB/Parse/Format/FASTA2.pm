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

# $Id: FASTA2.pm,v 1.18 2013/09/09 21:31:04 npb Exp $

###########################################################################
#
# Handles: fasta, tfastx  2.x
#
###########################################################################
package NPB::Parse::Format::FASTA2;

use NPB::Parse::Format::FASTA;
use NPB::Parse::Format::FASTA2::align;
use NPB::Parse::Format::FASTA2::lalign;

use strict;

use vars qw(
	    @ISA

	    @VERSIONS

	    $NULL

	    $ENTRY_START
	    $ENTRY_END

	    $HEADER_START  
	    $HEADER_END    
	                   
	    $RANK_START 
	    $RANK_END   
	                   
	    $TRAILER_START 
	    $TRAILER_END   
	                   
	    $MATCH_START     
	    $MATCH_END       
	                   
	    $SUM_START     
	    $SUM_END       
	                   
	    $ALN_START     
	    $ALN_END       
);

@ISA   = qw(NPB::Parse::Format::FASTA);

@VERSIONS = ( 
	     '2' => [
		     'FASTA',
		     'TFASTX',
		     'ALIGN',
		     'LALIGN',
		     'SSEARCH',
		    ],
	    );

$NULL  = '^\s*$';#for emacs';

$ENTRY_START   = '(?:'
    . '^\s*FASTA searches a protein or DNA sequence data bank'
    . '|'
    . '^\s*TFASTX translates and searches a DNA sequence data bank'
    . '|'
    . '^\s*\S+\s*[,:]\s+\d+\s+(?:aa|nt)'
    . '|'
    . 'SSEARCH searches a sequence database'
    . '|'
    . $NPB::Parse::Format::FASTA2::align::ALIGN_START
    . '|'
    . $NPB::Parse::Format::FASTA2::lalign::ALIGN_START
    . ')';
$ENTRY_END     = '(?:'
    . 'Library scan:'
    . '|'
    . $NPB::Parse::Format::FASTA2::align::ALIGN_END
    . '|'
    . $NPB::Parse::Format::FASTA2::lalign::ALIGN_END
    . ')';
    
$HEADER_START  = $ENTRY_START;
$HEADER_END    = '^The best scores are:'; 
               
$RANK_START    = $HEADER_END;
$RANK_END      = $NULL;
               
$TRAILER_START = $ENTRY_END;
$TRAILER_END   = $ENTRY_END;

$MATCH_START   = '^>*\S{7}.*\(\d+ (?:aa|nt)\)';
$MATCH_END     = "(?:$MATCH_START|$ENTRY_END)";

$SUM_START     = $MATCH_START;
$SUM_END       = $NULL;
       
$ALN_START     = '^(?:\s+\d+\s+|\s+$)';  #the ruler
$ALN_END       = $MATCH_END;

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package NPB::Parse::Format::FASTA2::HEADER;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA::HEADER);

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

    $self->{'query'}     = '';
    $self->{'queryfile'} = '';

    while (defined ($line = $text->next_line)) {
    
	if ($line =~ /^\s*(version\s+(\S+).*)/) {
	    $self->{'full_version'} = $1;
	    $self->{'version'}      = $2;
	    next;
	}

	if ($line =~ /^\s*>?(\S+).*\:\s+(\d+)\s+(?:aa|nt)/) {
	    if ($self->{'queryfile'} eq '') {
		$self->{'queryfile'} = $1;
		$self->{'queryfile'} =~ s/,$//;
	    } else {
		$self->test_args($line, $1);
		$self->{'query'} = NPB::Parse::Record::clean_identifier($1);
	    }
	    $self->{'length'} = $2;
	    next;
	}

	if ($line =~ /^(\d+)\s+residues\s+in\s+(\d+)\s+sequences/) {
	    $self->test_args($line, $1,$2);
	    (
	     $self->{'residues'},
	     $self->{'sequences'},
	    ) = ($1, $2);
	    next;
	} 

	#ignore any other text

    }

    if (! defined $self->{'full_version'} ) {
	#can't determine version: hardwire one!
	$self->{'full_version'} = 'looks like FASTA 2';
	$self->{'version'}      = '2';
    }

    $self;
}


###########################################################################
package NPB::Parse::Format::FASTA2::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA1::MATCH::ALN);


###########################################################################
1;
