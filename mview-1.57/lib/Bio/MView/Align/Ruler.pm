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

# $Id: Ruler.pm,v 1.12 2005/12/12 20:42:48 brown Exp $

###########################################################################
package Bio::MView::Align::Ruler;

use Bio::MView::Align;
use Bio::MView::Display;
use Bio::MView::Align::Row;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Align::Row);

sub new {
    my $type = shift;
    my ($length, $subtype) = (@_, 'ruler');
    my $self = {};
    bless $self, $type;

    $self->{'length'} = $length;
    $self->{'type'}   = $subtype;

    $self->reset_display;

    $self;
}

sub id     { $_[0] }
sub string { '' }
sub length { $_[0]->{'length'} }

sub reset_display { 
    $_[0]->{'display'} =
	{
	 'type'     => 'ruler',
	 'label0'   => '',
	 'label1'   => '',
	 'label2'   => '',
	 'label3'   => '',
	 'label4'   => '',
	 'label5'   => '',
	 'label6'   => '',
	 'range'    => [],
	 'number'   => 1,
	};
    $_[0];
}


###########################################################################
1;
