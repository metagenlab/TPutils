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

# $Id: lalign.pm,v 1.4 2013/09/09 21:31:05 npb Exp $

###########################################################################
# FASTA2 suite LALIGN pairwise method.
###########################################################################
package NPB::Parse::Format::FASTA2::lalign;

use NPB::Parse::Format::FASTA;
use strict;

use vars qw(
	    @ISA
	    
	    $ALIGN_START
	    $ALIGN_END
);

@ISA = qw(NPB::Parse::Format::FASTA);

#delimit full ALIGN or LALIGN entry
$ALIGN_START             = '^\s*(?:LALIGN finds|Comparison of:)';
$ALIGN_END               = $ALIGN_START;

#ALIGN record types
my $ALIGN_ALNSEP         = '^----------';
my $ALIGN_ALIGNMENT      = '^(?:\s+\d+\s+|\S+)';    #the ruler or any text
my $ALIGN_ALIGNMENTend   = "(?:$ALIGN_ALNSEP|$ALIGN_END)";
my $ALIGN_SUMMARY        = '^\s*\S+\s+identity';
my $ALIGN_SUMMARYend     = $ALIGN_ALIGNMENT;
my $ALIGN_MATCH          = $ALIGN_SUMMARY;
my $ALIGN_MATCHend       = $ALIGN_ALIGNMENTend;
my $ALIGN_HEADER         = $ALIGN_START;
my $ALIGN_HEADERend      = $ALIGN_MATCH;
my $ALIGN_Null           = '^\s*$';#'

#Parse one entry
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

	#HEADER lines
	if ($line =~ /$ALIGN_HEADER/o) {
	    $text->scan_until($ALIGN_HEADERend, 'HEADER');
	    next;
	}
	
	#consume data

	#MATCH lines		       	      
	if ($line =~ /$ALIGN_MATCH/o) {
	    $text->scan_until($ALIGN_MATCHend, 'MATCH');
	    next;			       	      
	}				       	      
	
	#skip LALIGN alignment delimiter
	next    if $line =~ /$ALIGN_ALNSEP/o;

	#blank line or empty record: ignore
	next    if $line =~ /$ALIGN_Null/o;

	#default
	$self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::FASTA2::lalign::HEADER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

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

    $self->{'program'} = '?';
    $self->{'version'} = '?';
    $self->{'moltype'} = '?';

    $self->{'id1'}     = '?';
    $self->{'desc1'}   = '';
    $self->{'length1'} = 0;

    $self->{'id2'}     = '?';
    $self->{'desc2'}   = '';
    $self->{'length2'} = 0;

    #consume Name lines
    while (defined ($line = $text->next_line)) { 
	#print $line;	

	#program information
	if ($line =~ /^LALIGN/x) {
	    $self->{'program'} = 'LALIGN';
	    next;
	}
	
	#guess program information for LALIGN missing stderr
	if ($line =~ /^\s*Comparison of:/) {
	    if ($self->{'program'} eq '?') {
		$self->{'program'} = 'LALIGN';
	    }
	    next;
	}

	if ($line =~ /^
	    \s*version\s+(\S+)
	    /xo) {

	    $self->test_args($line, $1);
	    (
	     $self->{'version'},
	    ) = ($1);
	    next;
	}

	#first/second sequence: id, description, length
	if ($line =~ /^
	    \s*\(.\)\s+\S+\s+(\S+)\s+	#id
	    \s*(.*)\s	                #description (empty?)
	    -\s*(\d+)\s+(aa|nt)             #length
	    /xo) {
	    
	    if ($self->{'id1'} eq '?') {
		$self->test_args($line, $1, $3, $4);
		(
		 $self->{'id1'},
		 $self->{'desc1'},
		 $self->{'length1'},
		 $self->{'moltype'},
		 ) = ($1, (defined $2?$2:''), $3, $4);
	    } else {
		$self->test_args($line, $1, $3);
		(
		 $self->{'id2'},
		 $self->{'desc2'},
		 $self->{'length2'},
		 ) = ($1, (defined $2?$2:''), $3);
	    }
	    next;
	}

	#ignore any other text
    }
    $self;
}

sub print {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    NPB::Parse::Record::print $self, $indent;
    printf "$x%20s -> %s\n",   'program', $self->{'program'};
    printf "$x%20s -> %s\n",   'version', $self->{'version'};
    printf "$x%20s -> %s\n",   'moltype', $self->{'moltype'};
    printf "$x%20s -> %s\n",   'id1',     $self->{'id1'};
    printf "$x%20s -> %s\n",   'desc1',   $self->{'desc1'};
    printf "$x%20s -> %s\n",   'length1', $self->{'length1'};
    printf "$x%20s -> %s\n",   'id2',     $self->{'id2'};
    printf "$x%20s -> %s\n",   'desc2',   $self->{'desc2'};
    printf "$x%20s -> %s\n",   'length2', $self->{'length2'};
}


###########################################################################
package NPB::Parse::Format::FASTA2::lalign::MATCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    my $type = shift;
    if (@_ < 2) {
	#at least two args, ($offset, $bytes are optional).
	NPB::Message::die($type, "new() invalid argument list (@_)");
    }
    my ($parent, $text, $offset, $bytes) = (@_, -1, -1);
    my ($self, $line, $record);

    $self = new NPB::Parse::Record($type, $parent, $text, $offset, $bytes);
    $text = new NPB::Parse::Record_Stream($self);

    while (defined ($line = $text->next_line)) {

	if ($line =~ /$ALIGN_SUMMARY/o) {
	    $text->scan_until($ALIGN_SUMMARYend, 'SUM');
	    next;
	}

	if ($line =~ /$ALIGN_ALIGNMENT/o) {
	    $text->scan_until($ALIGN_ALIGNMENTend, 'ALN');
	    next;
	}

	#blank line or empty record: ignore
        next    if $line =~ /$ALIGN_Null/o;

	#default
	$self->warn("unknown field: $line");
    }
    $self;
}

###########################################################################
package NPB::Parse::Format::FASTA2::lalign::MATCH::SUM;

use NPB::Parse::Regexps;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

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

    $self->{'identity'} = '?';
    $self->{'score'}    = '?';

    #consume Name lines
    while (defined ($line = $text->next_line)) { 
	#print $line;

	#LALIGN
	if ($line =~ /^
	    \s*($RX_Ureal)%\s+identity
	    .*
	    score:\s+($RX_Sint)
	    /xo) {

	    $self->{'identity'} = $1;
	    $self->{'score'}    = $2;
	    next;
	}

	#blank line or empty record: ignore
        next    if $line =~ /$ALIGN_Null/o;

	#default
	$self->warn("unknown field: $line");
    }
    $self;
}

sub print {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    NPB::Parse::Record::print $self, $indent;
    printf "$x%20s -> %s\n",   'identity', $self->{'identity'};
    printf "$x%20s -> %s\n",   'score',    $self->{'score'};
}

###########################################################################
package NPB::Parse::Format::FASTA2::lalign::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA2::MATCH::ALN);


###########################################################################
1;
