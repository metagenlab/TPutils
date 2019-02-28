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

# $Id: fastx.pm,v 1.2 2013/09/09 21:31:05 npb Exp $

###########################################################################
package NPB::Parse::Format::FASTA3::fastx;

use NPB::Parse::Format::FASTA3;
use NPB::Parse::Format::FASTA3::fasta;
use strict;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package NPB::Parse::Format::FASTA3::fastx::HEADER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3::fasta::HEADER);


###########################################################################
package NPB::Parse::Format::FASTA3::fastx::RANK;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3::fasta::RANK);


###########################################################################
package NPB::Parse::Format::FASTA3::fastx::TRAILER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3::fasta::TRAILER);


###########################################################################
package NPB::Parse::Format::FASTA3::fastx::MATCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3::fasta::MATCH);


###########################################################################
package NPB::Parse::Format::FASTA3::fastx::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3::fasta::MATCH::SUM);


###########################################################################
package NPB::Parse::Format::FASTA3::fastx::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3::MATCH::ALN);

# fast[xy]  dna x pro

sub query_base { return 3 }


###########################################################################
1;
