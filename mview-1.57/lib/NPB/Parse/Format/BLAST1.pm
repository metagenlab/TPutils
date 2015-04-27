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

# $Id: BLAST1.pm,v 1.12 2013/09/09 21:31:04 npb Exp $

###########################################################################
#
# Base classes for NCBI BLAST1, WashU BLAST2 families.
#
# Handles: BLAST 1.4.x, WashU 2.0x
#
# BLAST (pre NCBI version 2) parsing consists of 6 main record types:
#
#   HEADER        the header text
#   WARNING       optional warning messages
#   HISTOGRAM     the optional scores histogram
#   RANK          the list of ordered high scoring hits
#   MATCH         the set of fragments (HSPs) for a given hit
#   PARAMETERS    the trailer
#
# MATCH is further subdivided into:
#   SUM           the summary lines for each hit
#   ALN           each aligned fragment: score + alignment
#
###########################################################################
package NPB::Parse::Format::BLAST1;

use NPB::Parse::Format::BLAST;
use NPB::Parse::Regexps;

use strict;

use vars qw(@ISA

	    @VERSIONS

	    $NULL

	    $ENTRY_START
	    $ENTRY_END

            $WARNING_START
            $WARNING_END

            $WARNINGS_START
            $WARNINGS_END

            $PARAMETERS_START
            $PARAMETERS_END

            $MATCH_START
            $MATCH_END

            $RANK_START
            $RANK_MATCH
            $RANK_END

            $HISTOGRAM_START
            $HISTOGRAM_END

            $HEADER_START
            $HEADER_END

            $SCORE_START
            $SCORE_END
	   );

@ISA   = qw(NPB::Parse::Format::BLAST);

@VERSIONS = ( 
	     '1' => [
		     'BLASTP',
		     'BLASTN',
		     'BLASTX',
		     'TBLASTN',
		     'TBLASTX',
		    ],
	    );

$NULL  = '^\s*$';#for emacs';

$ENTRY_START      = '(?:'
    . '^BLASTP'
    . '|'
    . '^BLASTN'
    . '|'
    . '^BLASTX'
    . '|'
    . '^TBLASTN'
    . '|'
    . '^TBLASTX'
    . ')';
$ENTRY_END        = '^WARNINGS\s+ISSUED:';

$WARNING_START	  = '^(?:WARNING|NOTE):';
$WARNING_END	  = $NULL;

$WARNINGS_START	  = '^WARNINGS\s+ISSUED:';
$WARNINGS_END	  = $NULL;

$PARAMETERS_START = '^Parameters';
$PARAMETERS_END	  = "(?:$WARNINGS_START|$ENTRY_START)";

#$MATCH_START	  = '(?:^>|^[^\|:\s]+\|)';  #fails for some long query text
$MATCH_START	  = '^>';                   #fails for web pages
$MATCH_END	  = "(?:$MATCH_START|$WARNING_START|$PARAMETERS_START|$PARAMETERS_END)";

$RANK_START	  = '^\s+Smallest';
$RANK_MATCH	  = "(?:$NULL|^\\s+Smallest|^\\s+Sum|^\\s+High|^\\s+Reading|^\\s*Sequences|^[^>].*$RX_Uint\\s+$RX_Ureal\\s+$RX_Uint|$NPB::Parse::Format::BLAST::GCG_JUNK)";
$RANK_END	  = "(?:$WARNING_START|$PARAMETERS_START|$PARAMETERS_END)";

$HISTOGRAM_START  = '^\s+Observed Numbers';
$HISTOGRAM_END	  = "(?:$RANK_START|$RANK_END)";

$HEADER_START     = $ENTRY_START;
$HEADER_END       = "(?:$WARNING_START|$HISTOGRAM_START|$HISTOGRAM_END)";

$SCORE_START      = '^ Score';
$SCORE_END        = "(?:$SCORE_START|$MATCH_START|$PARAMETERS_START)";


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

	#Histogram lines
	if ($line =~ /$HISTOGRAM_START/o) {
	    $text->scan_until($HISTOGRAM_END, 'HISTOGRAM');
	    next;
	}

	#Rank lines: override $RANK_END definition
	if ($line =~ /$RANK_START/o) {       	      
	    $text->scan_while($RANK_MATCH, 'RANK');
	    next;
	}				       	      
	
	#Hit lines
	if ($line =~ /$MATCH_START/o) {
	    $text->scan_until($MATCH_END, 'MATCH');
	    next;			       	      
	}				       	      
	
	#WARNING lines
	if ($line =~ /$WARNING_START/o) {       	      
	    $text->scan_until($WARNING_END, 'WARNING');
	    next;			       	      
	}				       	      

	#Parameter lines
	if ($line =~ /$PARAMETERS_START/o) {       	      
	    $text->scan_until($PARAMETERS_END, 'PARAMETERS');
	    next;			       	      
	}				       	      
	
	#WARNINGS ISSUED line: ignore
	next    if $line =~ /$WARNINGS_START/o;

	#blank line or empty record: ignore
	next    if $line =~ /$NULL/o;

	#default
	$self->warn("unknown field: $line");
    }

    $self;#->examine;
}


###########################################################################
1;
