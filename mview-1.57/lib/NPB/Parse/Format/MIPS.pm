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

# $Id: MIPS.pm,v 1.8 2005/12/12 20:42:48 brown Exp $

###########################################################################
package NPB::Parse::Format::MIPS;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);


#delimit full MIPS entry
#my $MIPS_START          = '^\s*(?:\S+)?\s+MIPS:';
#my $MIPS_START          = '^(?:PileUp|\s*(?:\S+)?\s+MIPS:)';
my $MIPS_START          = '^>';
my $MIPS_END            = $MIPS_START;

#MIPS record types
my $MIPS_HEADER         = $MIPS_START;
my $MIPS_HEADERend      = '^L;';
my $MIPS_NAME           = $MIPS_HEADERend;
my $MIPS_NAMEend        = '^C;Alignment';
my $MIPS_ALIGNMENT      = '^\s*\d+';
my $MIPS_ALIGNMENTend   = $MIPS_START;
my $MIPS_Null           = '^\s*$';#'


#Consume one entry-worth of input on stream $fh associated with $file and
#return a new MIPS instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    my $fh   = $parent->{'fh'};
    my $text = $parent->{'text'};

    while (defined ($line = <$fh>)) {

	#start of entry
	if ($line =~ /$MIPS_START/o and $offset < 0) {
	    $offset = $fh->tell - length($line);
	    next;
	}

	#consume rest of stream
	last  if $line =~ /$MIPS_END/o;
    }
    return 0   if $offset < 0;

    $bytes = $fh->tell - $offset;

    new NPB::Parse::Format::MIPS(undef, $text, $offset, $bytes);
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
	if ($line =~ /$MIPS_HEADER/o) {
	    $text->scan_until($MIPS_HEADERend, 'HEADER');
	    next;
	}
	
	#consume data

	#NAME lines	       	      
	if ($line =~ /$MIPS_NAME/o) {    
	    $text->scan_until($MIPS_NAMEend, 'NAME');
	    next;			       	      
	}				       	      
	
	#ALIGNMENT lines		       	      
	if ($line =~ /$MIPS_ALIGNMENT/o) {       	      
	    $text->scan_until($MIPS_ALIGNMENTend, 'ALIGNMENT');
	    next;			       	      
	}				       	      
	
	#blank line or empty record: ignore
	next    if $line =~ /$MIPS_Null/o;

	#end of NAME section: ignore
	next    if $line =~ /$MIPS_NAMEend/o;
	
	#default
	$self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::MIPS::HEADER;

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

    $self->{'desc'} = '';

    #consume Name lines
    while (defined ($line = $text->next_line)) { 

	#> line
	if ($line =~ /^>[^;]+;(\S+)/o) {
	    $self->test_args($line, $1);
	    $self->{'ac'} = $1;
	    next;
	}
	
	#accumulate other lines
	$self->{'desc'} .= $line;
    }

    $self->warn("missing MIPS data\n")  unless exists $self->{'ac'};

    $self->{'desc'} = NPB::Parse::Record::strip_english_newlines($self->{'desc'});

    $self;
}

sub print {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    NPB::Parse::Record::print $self, $indent;
    printf "$x%20s -> %s\n",   'ac',     $self->{'ac'};
    printf "$x%20s -> '%s'\n", 'desc',   $self->{'desc'};
}


###########################################################################
package NPB::Parse::Format::MIPS::NAME;

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

    $self->{'seq'}   = {};
    $self->{'order'} = [];
    
    #consume Name lines
    while (defined ($line = $text->next_line)) {

	if ($line =~ /^L;(\S+)\s+(.*)/o) {
	    $self->test_args($line, $1,$2);
	    $self->{'seq'}->{$1} = NPB::Parse::Record::strip_english_newlines($2);
	    push @{$self->{'order'}}, $1;
	    next;
	} 
	
	next  if $line =~ /$MIPS_Null/;

	#default
	$self->warn("unknown field: $line");
    }
    $self;
}

sub print {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    local $_;
    NPB::Parse::Record::print $self, $indent;
    foreach (@{$self->{'order'}}) {
	printf "$x%20s -> %-15s %s=%s\n", 
	'seq',    $_,
	'desc',   $self->{'seq'}->{$_};
    }
}


###########################################################################
package NPB::Parse::Format::MIPS::ALIGNMENT;

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

    local $^W=0;
    local $_;
    
    $self->{'seq'} = {};

    while (defined ($line = $text->next_line)) {
    
	no strict;

	#start/end positions
	next  if $line =~ /^\s*\d+[^0-9]*\d+\s*$/o;

	#id/sequence
	if ($line =~ /^\s*(\S+)\s+([^0-9]+)\s+\d+$/o) {
	    $self->test_args($line, $1, $2);
	    $self->{'seq'}->{$1} .= $2;
	    next;
	} 

	#default: ignore all other line types (site and consensus data)
    }

    foreach (keys %{$self->{'seq'}}) {
	$self->{'seq'}->{$_} =~ s/ //g;
    }

    $self;
}

sub print {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    NPB::Parse::Record::print $self, $indent;
    local $_;
    foreach (sort keys %{$self->{'seq'}}) {
	printf "$x%20s -> %-15s =  %s\n", 'seq', $_, $self->{'seq'}->{$_};
    }
}


###########################################################################
1;
