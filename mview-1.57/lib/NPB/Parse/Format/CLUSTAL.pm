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

# $Id: CLUSTAL.pm,v 1.19 2013/12/01 18:56:16 npb Exp $

###########################################################################
package NPB::Parse::Format::CLUSTAL;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);


#delimit full CLUSTAL entry
my $CLUSTAL_START          = '^\s*CLUSTAL';
my $CLUSTAL_END            = $CLUSTAL_START;

#CLUSTAL record types
my $CLUSTAL_ALIGNMENT      = '^\s*\S+\s+\S+(?:\s+\d+)?\s*$';
my $CLUSTAL_ALIGNMENTend   = $CLUSTAL_START;
my $CLUSTAL_HEADER         = $CLUSTAL_START;
my $CLUSTAL_HEADERend      = "(?:$CLUSTAL_ALIGNMENT|$CLUSTAL_HEADER)";
my $CLUSTAL_Null           = '^\s*$';#'


#Consume one entry-worth of input on stream $fh associated with $file and
#return a new CLUSTAL instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    my $fh   = $parent->{'fh'};
    my $text = $parent->{'text'};

    while (defined ($line = <$fh>)) {

	#start of entry
 	if ($line =~ /$CLUSTAL_START/o and $offset < 0) {
	    $offset = $fh->tell - length($line);
	    next;
	}

	#consume rest of stream
        if ($line =~ /$CLUSTAL_END/o) {
	    last;
        }
    }
    return 0   if $offset < 0;

    $bytes = $fh->tell - $offset;

    new NPB::Parse::Format::CLUSTAL(undef, $text, $offset, $bytes);
}

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
	if ($line =~ /$CLUSTAL_HEADER/o) {
	    $text->scan_until($CLUSTAL_HEADERend, 'HEADER');
	    next;
	}
	
	#consume data

	#ALIGNMENT lines		       	      
	if ($line =~ /$CLUSTAL_ALIGNMENT/o) {
	    $text->scan_until($CLUSTAL_ALIGNMENTend, 'ALIGNMENT');
	    next;			       	      
	}				       	      
	
	#blank line or empty record: ignore
	next    if $line =~ /$CLUSTAL_Null/o;

	#default
	$self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::CLUSTAL::HEADER;

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

    $self->{'version'} = '';
    $self->{'major'}   = '';
    $self->{'minor'}   = '';

    #consume Name lines
    while (defined ($line = $text->next_line)) { 

	#first part of CLUSTAL line
	if ($line =~ /^
	    \s*
	    CLUSTAL
	    \s+
	    (([^\(\s]+)    #major version, eg., W
	    \s*
	    \((\S+)\))     #minor version, eg., 1.70
	    /xo) {

	    $self->test_args($line, $1, $2, $3);
	    (
	     $self->{'version'},
	     $self->{'major'},
	     $self->{'minor'},
	    ) = ($1, $2, $3);
	    
	}
	
	#first part of T-COFFEE line
	if ($line =~ /^
	    \s*
	    CLUSTAL\sFORMAT\sfor\s((T-COFFEE)
	    \s+
	    Version_(\S+)),     #version number, eg., 1.37
	    /xo) {

	    $self->test_args($line, $1, $2, $3);
	    (
	     $self->{'version'},
	     $self->{'major'},
	     $self->{'minor'},
	    ) = ($1, $2, $3);
	    
	}
	
	#ignore any other text
    }

    $self;
}

sub print {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    NPB::Parse::Record::print $self, $indent;
    printf "$x%20s -> %s\n",   'version', $self->{'version'};
    printf "$x%20s -> %s\n",   'major',   $self->{'major'};
    printf "$x%20s -> %s\n",   'minor',   $self->{'minor'};
}


###########################################################################
package NPB::Parse::Format::CLUSTAL::ALIGNMENT;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    my $type = shift;
    if (@_ < 2) {
	#at least two args, ($offset, $bytes are optional).
	NPB::Message::die($type, "new() invalid arguments (@_)");
    }
    my ($parent, $text, $offset, $bytes) = (@_, -1, -1);
    my ($self, $id, $line, $record);
    
    $self = new NPB::Parse::Record($type, $parent, $text, $offset, $bytes);
    $text = new NPB::Parse::Record_Stream($self);

    local $^W=0;
    local $_;
    
    $self->{'id'}    = [];
    $self->{'seq'}   = {};
    $self->{'match'} = '';

    my $off = 0;

    while (defined ($line = $text->next_line)) {
    
	no strict;

	chomp $line;

	#match symbols, but only if expected
	if ($off and $line !~ /[^*:. ]/) {
	    $line = substr($line, $off);
	    $self->{'match'} .= $line;
	    $off = 0;
	    next;
	}

	#id/sequence (/optional number)
	if ($line =~ /^\s*(\S+)\s+(\S+)(?:\s+\d+)?$/o) {
	    $self->test_args($line, $1, $2);
	    push @{$self->{'id'}}, $1    unless exists $self->{'seq'}->{$1};
	    $self->{'seq'}->{$1} .= $2;
	    $off = length($line) - length($2);
	    next;
	} 

	next    if $line =~ /$CLUSTAL_Null/o;

	#default
	$self->warn("unknown field: $line");
    }

    #line length check (ignore 'match' as this may be missing)
    if (defined $self->{'id'}->[0]) {
	$off = length $self->{'seq'}->{$self->{'id'}->[0]};
	foreach $id (keys %{$self->{'seq'}}) {
	    $line = $self->{'seq'}->{$id};
	    my $len = length $line;
	    #warn "$off, $len, $id\n";
	    if ($len != $off) {
		$self->die("length mismatch for '$id' (expect $off, saw $len):\noffending sequence: [$line]\n");
	    }
	}
    }

    $self;
}

sub print {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    NPB::Parse::Record::print $self, $indent;
    local $_;
    foreach (@{$self->{'id'}}) {
	printf "$x%20s -> %-15s %s\n", 'seq', $_, $self->{'seq'}->{$_};
    }
    printf "$x%20s -> %-15s %s\n", 'match', '', $self->{'match'};
}


###########################################################################
1;
