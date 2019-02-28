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

# $Id: Universal.pm,v 1.26 2008/05/20 15:31:46 npb Exp $

######################################################################
package Universal;

#useful general stuff
$::Date = `date`;
$::Prog = basename($0);

use strict;

#sub member {
#    my ($pattern, $list) = @_;
#
#    if (scalar(grep(/^$pattern$/, @$list)) != 1) {
#       $pattern = '\\' . $pattern;
#       if (scalar(grep(/^$pattern$/, @$list)) != 1) {
#           return 0;
#       }
#    }
#    return 1;
#}

sub member {
    my ($pattern, $list) = @_;
    my $i;
    foreach $i (@$list) {
        return 1    if $i eq $pattern;
    }
    return 0;
}

#dump sorted list of all instance variables, or supplied list
sub examine {
    my $self = shift;
    my @keys = @_ ? @_ : sort keys %$self;
    my $key;
    print "Class $self\n";
    foreach $key (@keys) {
        printf "%16s => %s\n", $key,
	    defined $self->{$key} ? $self->{$key} : '';
    }
    $self;
}

#shallow copy
sub copy {
    my $self = shift;
    my $copy = {};
    local $_;
    foreach (keys %$self) {
	#warn "$_ => $self->{$_}\n";
	if (defined $self->{$_}) {
	    $copy->{$_} = $self->{$_};
	} else {
	    $copy->{$_} = '';
	}
    }
    bless $copy, ref $self;
}

#deep copy
sub deep_copy {
    my $self = shift;
    my $copy = {};
    my $type;
    local $_;
    foreach (keys %$self) {
	#warn "$_ => $self->{$_}\n";
	if (defined $self->{$_}) {
	    if ($type = ref $self->{$_}) {
		if (UNIVERSAL::member($type, { qw(SCALAR ARRAY HASH CODE)})) {
		    $copy->{$_} = $self->{$_};
		} else {
		    $copy->{$_} = $copy->{$_}->deep_copy;
		}
	    }
	    $copy->{$_} = $self->{$_};
	} else {
	    $copy->{$_} = '';
	}
    }
    bless $copy, ref $self;
}

#warn with error string
sub warn {
    my $self = shift;
    chomp $_[$#_];
    if (ref($self)) {
	warn "Warning ", ref($self), '::', @_, "\n";
	return;
    }
    warn "Warning ", $self, '::', @_, "\n";
}

#exit with error string
sub die {
    my $self = shift;
    chomp $_[$#_];
    if (ref($self)) {
	die "Died ", ref($self), '::', @_, "\n";
    }
    die "Died ", $self, '::', @_, "\n";
}

#replacement for /bin/basename
sub basename {
    my($path, $ext) = (@_, "");
    ($path) = "/$path" =~ /.*\/(.+)$/;
    if ($path =~ /(.*)$ext$/) {
        return $1;
    }
    $path;
}

#arithmetic min() function
sub min {
    my ($a, $b) = @_;
    $a < $b ? $a : $b;
}


#arithmetic max() function
sub max {
    my ($a, $b) = @_;
    $a > $b ? $a : $b;
}

#Linux only?
sub vmstat {
    my ($s) = (@_, '');
    local ($_, *TMP);
    if (open(TMP, "cat /proc/$$/stat|")) {
	$_=<TMP>; my @ps = split /\s+/; close TMP;
	CORE::warn sprintf "VMEM=%8gk  $s\n", $ps[22] / 1024;
    } else {
	CORE::warn sprintf "VMEM=?  $s\n";
    }
}


######################################################################
1;
