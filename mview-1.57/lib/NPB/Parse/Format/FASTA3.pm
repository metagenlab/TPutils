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

# $Id: FASTA3.pm,v 1.16 2013/09/09 21:31:04 npb Exp $

###########################################################################
#
# Handles: fasta,fastx,fasty,tfasta,tfastx,tfasty,tfastxy   3.x
#
###########################################################################
package NPB::Parse::Format::FASTA3;

use NPB::Parse::Format::FASTA;
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
	     '3' => [
		     'FASTA',
		     'FASTX',
		     'FASTY',
		     'TFASTA',
		     'TFASTX',
		     'TFASTY',
		     'TFASTXY',
		    ],
	    );

$NULL  = '^\s*$';#for emacs';

$ENTRY_START   = '(?:^\s*'
    . 'FASTA +searches +a +protein +or +DNA +sequence +data +bank'
    . '|'
    . 'FAST[XY] +compares +a +DNA +sequence +to +a +protein +sequence +data +bank'
    . '|'
    . 'TFAST[AXY] +compares +a +protein +to +a +translated +DNA +data +bank'
    . '|'
    . 'TFASTXY +compares +a +protein +to +a +translated +DNA +data +bank'
    . '|'
    . '\S+\s*[,:]\s+\d+\s+(?:aa|nt)'
    . ')';
$ENTRY_END     = 'Function used was';

$HEADER_START  = $ENTRY_START;
$HEADER_END    = '^The\s+best\s+scores\s+are:';  #3.2t05 added a space!
               
$RANK_START    = $HEADER_END;
$RANK_END      = $NULL;
               
$TRAILER_START = '^\d+\s+residues\s+in\s+\d+\s+query';
$TRAILER_END   = $ENTRY_END;

#$MATCH_START   = '^>*\S{7}.*\(\d+ (?:aa|nt)\)';
$MATCH_START   = '^>+\S{7}';  #fasta -L flag makes multiple description lines
$MATCH_END     = "(?:$MATCH_START|$TRAILER_START|$ENTRY_END)";

$SUM_START     = $MATCH_START;
$SUM_END       = $NULL;
       
$ALN_START     = '^(?:\s+\d+\s+|\s+$)';  #the ruler
$ALN_END       = $MATCH_END;

#Parse one entry: generic for all FASTA3
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

    while (defined ($line = $text->next_line)) {

	#Header lines
	if ($line =~ /$HEADER_START/o) {
	    $text->scan_until($HEADER_END, 'HEADER');
	    next;
	}

	#Rank lines		       	      
	if ($line =~ /$RANK_START/o) {
	    $text->scan_until($RANK_END, 'RANK');
	    next;
	}				       	      
	
	#Hit lines		       	      
	if ($line =~ /$MATCH_START/o) {
	    $text->scan_until($MATCH_END, 'MATCH');
	    next;
	}

	#Trailer lines
	if ($line =~ /$TRAILER_START/o) {
	    $text->scan_until_inclusive($TRAILER_END, 'TRAILER');
	    next;
	}
	
	#end of FASTA job
	next    if $line =~ /$ENTRY_END/o;
	
	#blank line or empty record: ignore
	next    if $line =~ /$NULL/o;

	#default
	$self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::FASTA3::HEADER;

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

	#fasta < 3.4
	if ($line =~ /^\s*([^\s>]+)\s+(\d+)\s+(?:aa|nt)/o) {
	    $self->test_args($line, $1, $2);
	    (
	     $self->{'queryfile'},
	     $self->{'length'},
	    ) = ($1, $2);
	    $self->{'queryfile'} =~ s/,$//;
	    next;
	} 

	#fasta 3.4
	if ($line =~ /^\s*Query library\s+(\S+)/o) {
	    $self->test_args($line, $1);
	    $self->{'queryfile'} = $1;
	    $self->{'queryfile'} =~ s/,$//;
	    next;
	} 
	if ($line =~ /^\s*\d+>+(\S+).*-\s+(\d+)\s+(?:aa|nt)/) {
	    $self->test_args($line, $1, $2);
	    (
	     $self->{'query'},
	     $self->{'length'},
	    ) = (NPB::Parse::Record::clean_identifier($1),
		 $2);
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
	$self->{'full_version'} = 'looks like FASTA 3';
	$self->{'version'}      = '3';
    }

    $self;
}


###########################################################################
package NPB::Parse::Format::FASTA3::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA1::MATCH::ALN);

## The FASTA 3.x series gradually switched over to supplying the numeric
## sequence ranges for the alignment in the alignment summary section,
## allowing simpler parsing of this information. However, the application of
## this seems to have been inconsistently applied in the different programs,
## so the old FASTA 1 and 2 series parser that extracts these positions from
## the alignment rulers was used instead of the newer FASTA 3X parser. Since
## the FASTA 3.x series is now obsolete, this is moot.


###########################################################################
1;
