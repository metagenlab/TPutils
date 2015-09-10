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

# $Id: MULTAS.pm,v 1.9 2005/12/12 20:42:48 brown Exp $

###########################################################################
#
# MULTAS parsing consists of repeat records of:
#
# BLOCK
#    LIST
#    ALIGNMENT
# 
#
###########################################################################
package NPB::Parse::Format::MULTAS;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);


#delimit full MULTAS entry
my $MULTAS_START          = '^Block';
my $MULTAS_END            = undef;

my $MULTAS_Null           = '^\s*$';#'

#MULTAS record types
my $MULTAS_BLOCK          = $MULTAS_START;
my $MULTAS_BLOCKend       = "(?:$MULTAS_BLOCK|$MULTAS_Null)";
my $MULTAS_LIST           = '^\s*\d+\s+seqs';
#my $MULTAS_LISTmid        = "^(?:SEED|USER>>)";
my $MULTAS_LISTmid        = "^(?:SEED|USER)";
my $MULTAS_LISTend        = "(?:(?:>>){0}|$MULTAS_BLOCKend)";
my $MULTAS_ALIGNMENT      = '(?:>>){0}';
my $MULTAS_ALIGNMENTend   = $MULTAS_BLOCKend;


#Consume one entry-worth of input on stream $fh associated with $file and
#return a new MULTAS instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    my $fh   = $parent->{'fh'};
    my $text = $parent->{'text'};

    while (defined ($line = <$fh>)) {
	
	#start of entry
	if ($line =~ /$MULTAS_START/o and $offset < 0) {
            $offset = $fh->tell - length($line);
	    next;
	}

	#end of entry
	#if ($line =~ /$MULTAS_END/o) {
	#    last;
	#}
    }
    return 0   if $offset < 0;

    $bytes = $fh->tell - $offset;

    new NPB::Parse::Format::MULTAS(undef, $text, $offset, $bytes);
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

	#BLOCK lines	       	      
	if ($line =~ /$MULTAS_BLOCK/o) {
	    $text->scan_until($MULTAS_BLOCKend, 'BLOCK');
	    next;			       	      
	}				       	      
	
	#blank line or empty record: ignore
	if ($line =~ /$MULTAS_Null/o) {
	    next;
	}

	#terminal line: ignore
	#if ($line =~ /$MULTAS_END/o) {
	#    next;
	#}

	#default
	$self->warn("unknown field: $line");
    }

    $self->test_records(qw(BLOCK));

    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::MULTAS::BLOCK;

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

    while (defined ($line = $text->next_line)) {

	#BLOCK number
	if ($line =~ /$MULTAS_BLOCK\s+(\d+)/) {
	    $self->{'number'} = $1;
	    next;
	}

	#LIST lines		       	      
	if ($line =~ /$MULTAS_LIST/o) {       	      
	    $text->scan_while($MULTAS_LISTmid, 'LIST');
	    next;			       	      
	}				       	      

#	#LIST lines		       	      
#	if ($line =~ /$MULTAS_LIST/o) {       	      
#	    $text->scan_until($MULTAS_LISTend, 'LIST');
#	    next;			       	      
#	}				       	      

	#ALIGNMENT lines		       	      
	if ($line =~ /$MULTAS_ALIGNMENT/o) {       	      
	    $text->scan_until($MULTAS_ALIGNMENTend, 'ALIGNMENT');
	    next;			       	      
	}				       	      

	#blank line or empty record: ignore
	if ($line =~ /$MULTAS_Null/o) {
	    next;
	}

	#default
	$self->warn("unknown field: $line");
    }

    $self->test_records(qw(LIST ALIGNMENT));

    $self;#->examine;
}

sub print {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    NPB::Parse::Record::print $self, $indent;
    printf "$x%20s -> %d\n",   'number',     $self->{'number'};
}


###########################################################################
package NPB::Parse::Format::MULTAS::BLOCK::LIST;

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

    $self->{'count'} = 0;
    $self->{'hit'}   = [];

    #ranked search hits
    while (defined ($line = $text->next_line)) {

	if ($line =~ /\s*(\d+)\s+seqs/) {
	    $self->{'count'} = $1;
	    next;
	}

	if ($line =~ /\s*
	    ((?:SEED|USER))              #source of sequence
	    #\s*>>\s*
	    \s*>+\s*
	    ([^=]+)                      #identifier
	    \s*=\s*
	    (.*)                         #description
	    /xo) {

	    $self->test_args($line, $1, $2, $3);

	    push @{$self->{'hit'}},
	    {
	     'type'    => $1,
	     'id'      => $2,
	     'desc'    => $3,
	    };

	    #strip leading, trailing, internal white space
	    $self->{'hit'}->[$#{$self->{'hit'}}]->{'id'} =~ s/\s//g;

	    next;
	}
	
	#default
	$self->warn("unknown field: $line");
    }

    if ($self->{'count'} != @{$self->{'hit'}}) {
	$self->warn("stated ($self->{'count'}) and parsed (@{[scalar @{$self->{'hit'}}]}) sequence count mismatch");
    }

    $self;
}

sub print {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my ($hit, $field);
    NPB::Parse::Record::print $self, $indent;
    printf "$x%20s -> %d\n", 'count',        $self->{'count'};
    foreach $hit (@{$self->{'hit'}}) {
	foreach $field (sort keys %$hit) {
	    printf "$x%20s -> '%s'\n", $field,  $hit->{$field};
	}
    }
}


###########################################################################
package NPB::Parse::Format::MULTAS::BLOCK::ALIGNMENT;

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

    my (@tmp, $i);

    $self->{'seq'}    = [];
    $self->{'length'} = 0;
    $self->{'count'}  = 0;

    while (defined ($line = $text->next_line)) {

	chomp $line; @tmp = split(//, $line);
	
	for ($i=0; $i<@tmp; $i++) {
	    push @{$self->{'seq'}->[$i]}, $tmp[$i];
	}
    }

    $self->{'length'} = @{$self->{'seq'}->[0]};
    $self->{'count'}  = @{$self->{'seq'}};

    #check alignments all same length, then make string
    foreach $i (@{$self->{'seq'}}) {
	if (@$i != $self->{'length'}) {
	    $self->warn("alignment lengths differ");
	}
	$i = join('', @$i);
    }

    $self;
}

sub print {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my ($align, $i);
    NPB::Parse::Record::print $self, $indent;
    printf "$x%20s -> %d\n",    'count',  $self->{'count'};
    printf "$x%20s -> %d\n",    'length', $self->{'length'};
    for ($i=0; $i<@{$self->{'seq'}}; $i++) {
	printf("$x%20s[%d] -> '%s'\n", 'seq', $i+1, $self->{'seq'}->[$i]);
    }
}


###########################################################################
1;
