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

# $Id: tblastn.pm,v 1.8 2005/12/12 20:42:48 brown Exp $

###########################################################################
package NPB::Parse::Format::BLAST1::tblastn;

use NPB::Parse::Format::BLAST1::blastx;
use NPB::Parse::Format::BLAST1::tblastx;
use strict;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST1::blastx);


###########################################################################
package NPB::Parse::Format::BLAST1::tblastn::HEADER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST1::blastx::HEADER);


###########################################################################
package NPB::Parse::Format::BLAST1::tblastn::RANK;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST1::tblastx::RANK);


###########################################################################
package NPB::Parse::Format::BLAST1::tblastn::MATCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST1::blastx::MATCH);


###########################################################################
package NPB::Parse::Format::BLAST1::tblastn::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST1::blastx::MATCH::SUM);


###########################################################################
package NPB::Parse::Format::BLAST1::tblastn::MATCH::ALN;

use vars qw(@ISA);
use NPB::Parse::Regexps;

@ISA = qw(NPB::Parse::Format::BLAST::MATCH::ALN);

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

    #Score line
    $line = $text->next_line;

    if ($line =~ /^\s*
	Score\s*=\s*
	($RX_Uint)                           #score
	\s+
	\(($RX_Ureal)\s+bits\),              #bits
	\s+
	Expect\s*=\s*
	($RX_Ureal),                         #expectation
	\s+
	(?:Sum\sP\((\d+)\)|P)\s*=\s*         #number of frags
	($RX_Ureal)                          #p-value
	/xo) {
	
	$self->test_args($line, $1, $2, $3, $5);

	(
	 $self->{'score'},
	 $self->{'bits'},
	 $self->{'expect'},
	 $self->{'n'},                       #substitute 1 unless $4
	 $self->{'p'},
	) = ($1, $2, $3, defined $4?$4:1, $5);
    }
    else {
	$self->warn("expecting 'Score' line: $line");
    }
    
    #Identities line
    $line = $text->next_line;

    if ($line =~ /^\s*
	Identities\s*=\s*
	(\d+\/\d+)                           #identities fraction
	\s+
	\((\d+)%\),                          #identities percentage
	\s+
	Positives\s*=\s*
	(\d+\/\d+)                           #positives fraction
	\s+
	\((\d+)%\)                           #positives percentage
	,\s+Frame\s*=\s*
	([+-])                               #frame sign
	(\d+)                                #frame number
	/xo) {
	
	$self->test_args($line, $1, $2, $3, $4, $5, $6);

	(
	 $self->{'id_fraction'},
	 $self->{'id_percent'},
	 $self->{'pos_fraction'},
	 $self->{'pos_percent'},
	 $self->{'sbjct_frame'},
	) = ($1, $2, $3, $4, $5 . $6);
	
	#record query orientation in MATCH list
	push @{$parent->{'orient'}->{$self->{'sbjct_frame'}}}, $self;
	
    } else {
	$self->warn("expecting 'Identities' line: $line");
    }

    $self->parse_alignment($text);

    $self;
}

sub print {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    NPB::Parse::Record::print $self, $indent;
    printf "$x%20s -> %s\n",   'score',          $self->{'score'};
    printf "$x%20s -> %s\n",   'bits',           $self->{'bits'};
    printf "$x%20s -> %s\n",   'expect',         $self->{'expect'};
    printf "$x%20s -> %s\n",   'p',              $self->{'p'};
    printf "$x%20s -> %s\n",   'n',              $self->{'n'};
    printf "$x%20s -> %s\n",   'id_fraction',    $self->{'id_fraction'};
    printf "$x%20s -> %s\n",   'id_percent',     $self->{'id_percent'};
    printf "$x%20s -> %s\n",   'pos_fraction',   $self->{'pos_fraction'};
    printf "$x%20s -> %s\n",   'pos_percent',    $self->{'pos_percent'};
    printf "$x%20s -> %s\n",   'sbjct_frame',    $self->{'sbjct_frame'};
    $self->SUPER::print($indent);
}


###########################################################################
package NPB::Parse::Format::BLAST1::tblastn::WARNING;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST1::blastx::WARNING);


###########################################################################
package NPB::Parse::Format::BLAST1::tblastn::HISTOGRAM;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST1::blastx::HISTOGRAM);


###########################################################################
package NPB::Parse::Format::BLAST1::tblastn::PARAMETERS;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST1::blastx::PARAMETERS);


###########################################################################
1;
