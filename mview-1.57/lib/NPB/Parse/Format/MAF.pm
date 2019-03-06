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

# $Id: MAF.pm,v 1.1 2013/12/01 18:56:16 npb Exp $

###########################################################################
package NPB::Parse::Format::MAF;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);


#delimit full MAF entry
my $MAF_START          = '^##maf\s';
my $MAF_END            = $MAF_START;

#MAF record types
my $MAF_Null           = '^\s*$';#'
my $MAF_BLOCK          = '^a\s';
my $MAF_BLOCKend       = $MAF_Null;
my $MAF_HEADER         = $MAF_START;
my $MAF_HEADERend      = "(?:$MAF_BLOCK|$MAF_HEADER)";


#Consume one entry-worth of input on stream $fh associated with $file and
#return a new MAF instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    my $fh   = $parent->{'fh'};
    my $text = $parent->{'text'};

    while (defined ($line = <$fh>)) {

	#start of entry
 	if ($line =~ /$MAF_START/o and $offset < 0) {
	    $offset = $fh->tell - length($line);
	    next;
	}

	#consume rest of stream
        if ($line =~ /$MAF_END/o) {
	    last;
        }
    }
    return 0   if $offset < 0;

    $bytes = $fh->tell - $offset;

    new NPB::Parse::Format::MAF(undef, $text, $offset, $bytes);
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
	if ($line =~ /$MAF_HEADER/o) {
	    $text->scan_until($MAF_HEADERend, 'HEADER');
	    next;
	}
	
	#consume data

	#BLOCK lines
	if ($line =~ /$MAF_BLOCK/o) {
	    $text->scan_until($MAF_BLOCKend, 'BLOCK');
	    next;
	}

	#blank line or empty record: ignore
	next    if $line =~ /$MAF_Null/o;

	#default
	$self->warn("unknown field: $line");
    }

    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::MAF::HEADER;

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

    #consume Name lines
    while (defined ($line = $text->next_line)) {

	#first part of MAF line
	if ($line =~ /^
            ##maf
	    \s+
	    version=(\S+)
	    /xo) {

	    $self->test_args($line, $1);
	    (
	     $self->{'version'},
	    ) = ($1);

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
}


###########################################################################
package NPB::Parse::Format::MAF::BLOCK;

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

    $self->{'row'} = [];

    my $off = 0;

    while (defined ($line = $text->next_line)) {

	no strict;

	chomp $line;

	#a score=xxxx.yyyy
	if ($line =~ /^a\sscore=(\S+)/) {
	    $self->test_args($line, $1);
	    $self->{'score'} = $1;
	    next;
	}

	#s src start size strand srcsize sequence
	if ($line =~ /^s\s(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s+(\d+)\s+(\S+)$/) {
	    $self->test_args($line, $1, $2, $3, $4, $5, $6);
	    push @{$self->{'row'}}, {
		'id'      => $1,
		'start'   => $2,
		'size'    => $3,
		'strand'  => $4,
		'srcsize' => $5,
		'seq'     => $6,
	    };
	    $off = length($line) - length($6);
	    next;
	}

	next    if $line =~ /^[aie]\s/o;
	next    if $line =~ /$MAF_Null/o;

	#default
	$self->warn("unknown field: $line");
    }

    #line length check
    if (defined $self->{'row'}->[0]) {
	$off = length $self->{'row'}->[0]->{'seq'};
	foreach my $row (@{$self->{'row'}}) {
	    $line = $row->{'seq'};
	    my $len = length $line;
	    #warn "$off, $len, $row->{'id'}\n";
	    if ($len != $off) {
		$self->die("length mismatch for '$row->{'id'}' (expect $off, saw $len):\noffending sequence: [$line]\n");
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
    printf "$x%20s -> %s\n", 'score', $self->{'score'};
    foreach my $row (@{$self->{'row'}}) {
	printf "$x%20s -> %s\n", 'id',      $row->{'id'};
	printf "$x%20s -> %s\n", 'start',   $row->{'start'};
	printf "$x%20s -> %s\n", 'size',    $row->{'size'};
	printf "$x%20s -> %s\n", 'strand',  $row->{'strand'};
	printf "$x%20s -> %s\n", 'srcsize', $row->{'srcsize'};
	printf "$x%20s -> %s\n", 'seq',     $row->{'seq'};
    }
}


###########################################################################
1;
