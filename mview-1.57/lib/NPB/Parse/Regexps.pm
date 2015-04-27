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

# $Id: Regexps.pm,v 1.7 2005/12/12 20:42:48 brown Exp $

###########################################################################
# regexps for string matching numerical types
###########################################################################
package NPB::Parse::Regexps;

use Exporter;

@ISA = qw(Exporter);

@EXPORT = 
    qw(
       $RX_Uint
       $RX_Sint
       $RX_Ureal 
       $RX_Sreal
      );


#unsigned integer
$RX_Uint   = '\+?\d+';

#signed integer
$RX_Sint   = '[+-]?\d+';

#unsigned real
$RX_Ureal = '\+?(?:\d+\.\d+|\d+\.|\d+|\.\d+)?(?:[eE][+-]?\d+)?';

#signed real
$RX_Sreal = '[+-]?(?:\d+\.\d+|\d+\.|\d+|\.\d+)?(?:[eE][+-]?\d+)?';


###########################################################################
1;
