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

# $Id: Pearson.pm,v 1.15 2005/12/12 20:42:48 brown Exp $

###########################################################################
package NPB::Parse::Format::Pearson;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);


#Pearson record types
my $Pearson_Null     = '^\s*$';#'
my $Pearson_SEQ      = '^\s*>';
my $Pearson_SEQend   = "(?:$Pearson_SEQ|$Pearson_Null)";


#Consume one entry-worth of input on stream $fh associated with $file and
#return a new Slurp instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    my $fh   = $parent->{'fh'};
    my $text = $parent->{'text'};

    while (defined ($line = <$fh>)) {
	
	#start of entry
	if ($offset < 0) {
            $offset = $fh->tell - length($line);
	    next;
	}

    }
    return 0   if $offset < 0;

    $bytes = $fh->tell - $offset;

    new NPB::Parse::Format::Pearson(undef, $text, $offset, $bytes);
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
	
	#SEQ lines		       	      
	if ($line =~ /$Pearson_SEQ/o) {
	    $text->scan_until($Pearson_SEQend, 'SEQ');
	    next;			       	      
	}				       	      
	
	#blank line or empty record: ignore
	if ($line =~ /$Pearson_Null/o) {
	    next;
	}

	#default
	$self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::Pearson::SEQ;

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

    $self->{'id'}    = '';
    $self->{'desc'}  = '';
    $self->{'seq'}   = '';
    
    while (defined ($line = $text->next_line(1))) {

	#read header line
	if ($line =~ /^\s*>\s*(\S+)\s*(.*)?/o) {
	    $self->test_args($line, $1);
	    (
	     $self->{'id'}, 
	     $self->{'desc'},
	    ) = ($1, "$2");
	    #2015-01-19, GeneDoc puts a '.' in after the identifier
	    $self->{'desc'} = ''  if $self->{'desc'} =~ /\s*\.\s*/;
	    next;
	} 

	#read sequence lines upto asterisk, if present
	if ($line =~ /([^\*]+)/) {
	    $self->{'seq'} .= $1;
	    next;
	}

	#ignore lone asterisk
	last    if $line =~ /\*/;

	#default
	$self->warn("unknown field: $line");
    }

    #strip internal whitespace from sequence
    $self->{'seq'} =~ s/\s//g;

    $self;
}

sub print {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    NPB::Parse::Record::print $self, $indent;
    printf "$x%20s -> %s\n",   'id',        $self->{'id'};
    printf "$x%20s -> '%s'\n", 'desc',      $self->{'desc'};
    printf "$x%20s -> %s\n",   'seq',       $self->{'seq'};
}


###########################################################################
1;
