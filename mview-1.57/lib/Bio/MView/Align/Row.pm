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

# $Id: Row.pm,v 1.23 2005/12/12 20:42:48 brown Exp $

###########################################################################
package Bio::MView::Align::Row;

use Bio::MView::Align;
use Bio::MView::Display;
use Bio::MView::Align::Ruler;
use Bio::MView::Align::Sequence;
use Bio::MView::Align::Identity;
use Bio::MView::Align::Consensus;
use strict;

use vars qw($Colour_Black $Colour_White $Colour_DarkGray $Colour_LightGray
	    $Colour_Cream $Colour_Comment);

$Colour_Black   	    	= '#000000';
$Colour_White   	    	= '#FFFFFF';
$Colour_DarkGray    	    	= '#666666';
$Colour_LightGray   	    	= '#999999';
$Colour_Cream            	= '#FFFFCC';
$Colour_Comment                 = '#aa6666';  #brown

#sub DESTROY { warn "DESTROY $_[0]\n" }

sub print {
    my $self = shift;
    local $_;
    foreach (sort keys %$self) {
	printf "%15s => %s\n", $_, $self->{$_};
    }
    $self;
}

sub set_display {
    my $self = shift;
    my ($key, $val, @tmp, %tmp);
    #$self->print;
    while ($key = shift @_) {

	$val = shift @_;

	#have to copy referenced data in case caller iterates over
	#many instances of self and passes the same data to each!
	if (ref $val eq 'HASH') {
	    %tmp = %$val;
	    $val = \%tmp;
	} elsif (ref $val eq 'ARRAY') {
	    @tmp = @$val;
	    $val = \@tmp;
	}
	#warn "($key) $self->{'display'}->{$key} --> $val\n";
	$self->{'display'}->{$key} = $val;
    }

    $self;
}

sub get_display { $_[0]->{'display'} }

sub color_special {}
sub color_by_type {}
sub color_by_identity {}
sub color_by_consensus_sequence {}
sub color_by_consensus_group {}


###########################################################################
1;
