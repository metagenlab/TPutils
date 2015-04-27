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

# $Id: Plain.pm,v 1.9 2005/12/12 20:42:48 brown Exp $

######################################################################
#
# The simplest input alignment is a column of id's and a column of aligned
# sequences of identical length with the entire alignment in one block.
#
###########################################################################
package NPB::Parse::Format::Plain;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);


my $Plain_Null           = '^\s*$';#'
my $Plain_Comment        = '^\s*\#';

#delimit full Plain entry
my $Plain_START          = '^\s*\S+\s+\S';
my $Plain_END            = $Plain_Null;

#Plain record types
my $Plain_ALIGNMENT      = '^\s*\S+\s+\S';
my $Plain_ALIGNMENTend   = $Plain_END;


#Consume one entry-worth of input on stream $fh associated with $file and
#return a new Plain instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    my $fh   = $parent->{'fh'};
    my $text = $parent->{'text'};

    while (defined ($line = <$fh>)) {

	#start of entry
 	if ($line =~ /$Plain_START/o and $offset < 0) {
	    $offset = $fh->tell - length($line);
	    next;
	}

	#end of entry
        if ($line =~ /$Plain_END/o and ! $offset < 0) {
	    last;
        }
    }
    return 0   if $offset < 0;

    $bytes = $fh->tell - $offset;

    new NPB::Parse::Format::Plain(undef, $text, $offset, $bytes);
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

	#consume data

	#comment: ignore
	next    if $line =~ /$Plain_Comment/;

	#ALIGNMENT lines		       	      
	if ($line =~ /$Plain_ALIGNMENT/o) {
	    $text->scan_until($Plain_ALIGNMENTend, 'ALIGNMENT');
	    next;			       	      
	}				       	      
	
	#blank line or empty record: ignore
	next    if $line =~ /$Plain_Null/o;

	#default
	$self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::Plain::ALIGNMENT;

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
    
    $self->{'id'}    = [];
    $self->{'seq'}   = {};
    $self->{'match'} = '';

    my $off = 0;

    while (defined ($line = $text->next_line)) {

	chomp $line;

	#id/sequence
	if ($line =~ /^\s*(\S+)\s+(\S+)\s*$/o) {
	    $self->test_args($line, $1, $2);
	    push @{$self->{'id'}}, $1    unless exists $self->{'seq'}->{$1};
	    $self->{'seq'}->{$1} .= $2;
	    $off = length($line) - length($2);
	    next;
	} 

	#default
	$self->warn("unknown field: $line");
    }

    #line length check (ignore 'match' as this may be missing)
    if (defined $self->{'id'}->[0]) {
	$off = length $self->{'seq'}->{$self->{'id'}->[0]};
	foreach $line (keys %{$self->{'seq'}}) {
	    $line = $self->{'seq'}->{$line};
	    if (length $line != $off) {
		$self->die("unequal line lengths (expect $off, got @{[length $line]})\n");
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
