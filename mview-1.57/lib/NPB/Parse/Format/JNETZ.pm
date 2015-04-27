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

# $Id: JNETZ.pm,v 1.8 2005/12/12 20:42:48 brown Exp $

###########################################################################
#
# JNET -z format parsing consists of:
#   ALIGNMENT
#
###########################################################################
package NPB::Parse::Format::JNETZ;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);

#delimit full JNETZ entry
my $JNETZ_START          = '^START PRED';
my $JNETZ_END            = '^END PRED';
my $JNETZ_Null           = '^\s*$';#'

my $JNETZ_ALIGNMENT      = $JNETZ_START;
my $JNETZ_ALIGNMENTend   = $JNETZ_END;


#Consume one entry-worth of input on stream $fh associated with $file and
#return a new JNETZ instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    my $fh   = $parent->{'fh'};
    my $text = $parent->{'text'};

    while (defined ($line = <$fh>)) {
	
	#start of entry
	if ($line =~ /$JNETZ_START/o and $offset < 0) {
            $offset = $fh->tell - length($line);
	    next;
	}

	#end of entry
	if ($line =~ /$JNETZ_END/o) {
	    last;
	}
    }
    return 0   if $offset < 0;

    $bytes = $fh->tell - $offset;

    new NPB::Parse::Format::JNETZ(undef, $text, $offset, $bytes);
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

	#ALIGNMENT lines	       	      
	if ($line =~ /$JNETZ_ALIGNMENT/o) {
	    $text->scan_until($JNETZ_ALIGNMENTend, 'ALIGNMENT');
	    next;			       	      
	}				       	      
	
	#blank line or empty record: ignore
	next  if $line =~ /$JNETZ_Null/o;

	#terminal line: ignore
	next  if $line =~ /$JNETZ_END/o;

	#default
	$self->warn("unknown field: $line");
    }

    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::JNETZ::ALIGNMENT;

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

    #initialise these:
    $self->{'query'} = [];
    $self->{'final'} = [];
    $self->{'conf'}  = [];
    $self->{'align'} = [];

    while (defined ($line = $text->next_line)) {

	if ($line =~
	    /^\s*(\S)
	    \s+(\S)
	    \s+\|
	    \s+(\S)
	    \s+.
	    \s+(\d)
	    /xo
	   ) {
	    push @{ $self->{'query'}}, $1;
	    push @{ $self->{'final'}}, $2;
	    push @{ $self->{'align'}}, $3;
	    push @{ $self->{'conf'}},  $4;
	    next;
	}

	#header line: ignore
	next  if $line =~ /$JNETZ_ALIGNMENT/;

	#blank line or empty record: ignore
	next  if $line =~ /$JNETZ_Null/o;

	#default
	$self->warn("unknown field: $line");
    }

    $self;#->examine;
}

sub print {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    NPB::Parse::Record::print $self, $indent;
    printf "$x%20s -> %s\n",   'query', scalar $self->get_query;
    printf "$x%20s -> %s\n",   'final', scalar $self->get_final;
    printf "$x%20s -> %s\n",   'conf',  scalar $self->get_conf;
    printf "$x%20s -> %s\n",   'align', scalar $self->get_align;
    $self;
}

sub get_query {
    my $self = shift;
    my @a = @{$self->{'query'}};
    return \@a    if wantarray;
    return join("", @a);
}

sub get_final {
    my $self = shift;
    my @a = @{$self->{'final'}};
    return \@a    if wantarray;
    return join("", @a);
}

sub get_conf {
    my $self = shift;
    my @a = @{$self->{'conf'}};
    return \@a    if wantarray;
    return join("", @a);
}

sub get_align {
    my $self = shift;
    my @a = @{$self->{'align'}};
    return \@a    if wantarray;
    return join("", @a);
}


###########################################################################
1;
