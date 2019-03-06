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

# $Id: Sequence.pm,v 1.24 2013/12/01 18:56:16 npb Exp $

###########################################################################
package Bio::MView::Align::Sequence;

use Bio::MView::Align;
use Bio::MView::Display;
use Bio::MView::Align::Row;
use strict;

use vars qw(@ISA
	    $Default_PRO_Colormap $Default_DNA_Colormap
            $Default_FIND_Colormap $Default_Colormap
	    %Template);

@ISA = qw(Bio::MView::Align::Row);

$Default_PRO_Colormap  = 'P1';    #default protein colormap name
$Default_DNA_Colormap  = 'D1';    #default nucleotide colormap name
$Default_FIND_Colormap = 'FIND';  #default find pattern colormap name
$Default_Colormap = $Default_PRO_Colormap;

%Template = 
    (
     'id'      => undef,     #identifier
     'from'    => undef,     #start number of sequence
     'type'    => undef,     #information about own subtype
     'string'  => undef,     #alignment string
     'display' => undef,     #hash of display parameters
    );

sub new {
    my $type = shift;
    #warn "${type}::new() (@_)\n";
    if (@_ < 2) {
	die "${type}::new() missing arguments\n";
    }
    my ($id, $string, $subtype) = (@_, 'sequence');

    my $self = { %Template };

    $self->{'id'}     = $id;
    $self->{'from'}   = $string->lo;
    $self->{'type'}   = $subtype;
    $self->{'string'} = $string;

    bless $self, $type;

    $self->reset_display;

    $self;
}

#sub DESTROY { warn "DESTROY $_[0]\n" }

sub print {
    sub _format {
	my ($self, $k, $v) = @_;
	$v = 'undef' unless defined $v;
	$v = "'$v'" if $v =~ /^\s*$/;
	return sprintf("  %-15s => %s\n", $k, $v)
    }
    my $self = shift;
    warn "$self\n";
    map { warn $self->_format($_, $self->{$_}) } sort keys %{$self};
    $self;
}

sub id     { $_[0]->{'id'} }
sub string { $_[0]->{'string'}->string }
sub array  { split //, $_[0]->{'string'}->string }
sub from   { $_[0]->{'from'} }
sub length { $_[0]->{'string'}->length }
sub ungapped_length { $_[0]->{'string'}->ungapped_length }

sub reset_display {
    $_[0]->{'display'} =
	{
	 'type'     => 'sequence',
	 'label1'   => $_[0]->{'id'},
	 'sequence' => $_[0]->{'string'},
	 'range'    => [],
	};
    $_[0];
}

sub get_color {
    my ($self, $c, $map) = @_;
    my ($index, $color, $trans);

    #warn "get_color: $c, $map";

    #set transparent(T)/solid(S)
    if (exists $Bio::MView::Align::Colormaps->{$map}->{$c}) {
	$trans = $Bio::MView::Align::Colormaps->{$map}->{$c}->[1];
	$index = $Bio::MView::Align::Colormaps->{$map}->{$c}->[0];
	$color = $Bio::MView::Align::Palette->[1]->[$index];

	#warn "$c $map\{$c} [$index] [$color] [$trans]\n";
	
	return ($color, "$trans$index");
    }

    #wildcard colour
    if (exists $Bio::MView::Align::Colormaps->{$map}->{'*'}) {

	$trans = $Bio::MView::Align::Colormaps->{$map}->{'*'}->[1];
	$index = $Bio::MView::Align::Colormaps->{$map}->{'*'}->[0];
	$color = $Bio::MView::Align::Palette->[1]->[$index];

	#warn "$c $map\{'*'} [$index] [$color] [$trans]\n";

	return ($color, "$trans$index");
    }

    #preset colour name in $map, used for string searches
    #where all matches should be same colour
    if (exists $Bio::MView::Align::Palette->[0]->{$map}) {

	$trans = 'S';
	$index = $Bio::MView::Align::Palette->[0]->{$map};
        $color = $Bio::MView::Align::Palette->[1]->[$index];

	#warn "$c $map\{$c} [$index] [$color] [$trans]\n";

	return ($color, "$trans$index");
    }

    return 0;    #no match
}

sub color_special {
    my $self = shift;
    my %par = @_;

    $par{'css1'}     = 0
	unless defined $par{'css1'};
    $par{'symcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'symcolor'};
    $par{'gapcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'gapcolor'};
    $par{'colormap'} = $Default_Colormap
	unless defined $par{'colormap'};

    #locate a 'special' colormap'
    my ($size, $map) = (0);
    foreach $map (keys %$Bio::MView::Align::Colormaps) {
	if ($self->{'id'} =~ /$map/i) {
	    if (length($&) > $size) {
		$par{'colormap'} = $map;
		$size = length($&);
	    }
	}
    }
    return unless $size;

    #warn $par{'colormap'};

    my ($color, $end, $i, $c, @tmp) = ($self->{'display'}->{'range'});
    
    push @$color, 1, $self->length, 'color' => $par{'symcolor'};

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

	$c = $self->{'string'}->raw($i);
	
	#warn "$self->{'id'}  [$i]= $c\n";

	#white space: no color
	next    if $self->{'string'}->is_space($c);

	#gap: gapcolour
	if ($self->{'string'}->is_non_sequence($c)) {
	    push @$color, $i, 'color' => $par{'gapcolor'};
	    next;
	}
	
	#use symbol color/wildcard colour
	@tmp = $self->get_color($c, $par{'colormap'});

	if (@tmp) {
	    if ($par{'css1'}) {
		push @$color, $i, 'class' => $tmp[1];
	    } else {
		push @$color, $i, 'color' => $tmp[0];
	    }
	} else {
	    push @$color, $i, 'color' => $par{'symcolor'};
	}
    }
    
    $self->{'display'}->{'paint'}  = 1;
    $self;
}

sub color_by_find_block {
    my $self = shift;
    my %par = @_;

    return unless $self->{'type'} eq 'sequence';

    $par{'css1'}     = 0
	unless defined $par{'css1'};
    $par{'symcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'symcolor'};
    $par{'gapcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'gapcolor'};
    $par{'colormap'} = $Default_Colormap
	unless defined $par{'colormap'};
    $par{'find'}     = ''
        unless defined $par{'find'};

    my ($color, $end, $i, $c, @tmp) = ($self->{'display'}->{'range'});
    
    push @$color, 1, $self->length, 'color' => $par{'symcolor'};

    my %mindex = ();
    foreach my $block (@{$self->{string}->findall($par{patterns},
						  $par{mapsize})}) {
	$mindex{$block->[1]} = $block->[0];
    }

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {
        
	$c = $self->{'string'}->raw($i);
	
	#warn "[$i]= $c\n";

	#white space: no color
	next    if $self->{'string'}->is_space($c);
        
	#gap: gapcolour
	if ($self->{'string'}->is_non_sequence($c)) {
	    push @$color, $i, 'color' => $par{'gapcolor'};
	    next;
	}
	
        if (exists $mindex{$i}) {
            #use symbol color/wildcard colour
            @tmp = $self->get_color($mindex{$i}, $par{'colormap'});
        } else {
            @tmp = ();
        }
	
	if (@tmp) {
	    if ($par{'css1'}) {
		push @$color, $i, 'class' => $tmp[1];
	    } else {
		push @$color, $i, 'color' => $tmp[0];
	    }
	} else {
	    push @$color, $i, 'color' => $par{'symcolor'};
	}
    }
    
    $self->{'display'}->{'paint'} = 1;
    $self;
}

sub color_by_type {
    my $self = shift;
    my %par = @_;

    return unless $self->{'type'} eq 'sequence';

    $par{'css1'}     = 0
	unless defined $par{'css1'};
    $par{'symcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'symcolor'};
    $par{'gapcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'gapcolor'};
    $par{'colormap'} = $Default_Colormap
	unless defined $par{'colormap'};

    my ($color, $end, $i, $c, @tmp) = ($self->{'display'}->{'range'});
    
    push @$color, 1, $self->length, 'color' => $par{'symcolor'};

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

	$c = $self->{'string'}->raw($i);
	
	#warn "[$i]= $c\n";

	#white space: no color
	next    if $self->{'string'}->is_space($c);

	#gap: gapcolour
	if ($self->{'string'}->is_non_sequence($c)) {
	    push @$color, $i, 'color' => $par{'gapcolor'};
	    next;
	}
	
	#use symbol color/wildcard colour
	@tmp = $self->get_color($c, $par{'colormap'});
	
	if (@tmp) {
	    if ($par{'css1'}) {
		push @$color, $i, 'class' => $tmp[1];
	    } else {
		push @$color, $i, 'color' => $tmp[0];
	    }
	} else {
	    push @$color, $i, 'color' => $par{'symcolor'};
	}
    }
    
    $self->{'display'}->{'paint'}  = 1;
    $self;
}

sub color_by_identity {
    my ($self, $othr) = (shift, shift);
    my %par = @_;

    return unless $self->{'type'} eq 'sequence';
    return unless $othr;

    die "${self}::color_by_identity() length mismatch\n"
	unless $self->length == $othr->length;

    $par{'css1'}     = 0
	unless defined $par{'css1'};
    $par{'symcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'symcolor'};
    $par{'gapcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'gapcolor'};
    $par{'colormap'} = $Default_Colormap
	unless defined $par{'colormap'};

    my ($color, $end, $i, $c1, $c2, @tmp) = ($self->{'display'}->{'range'});

    push @$color, 1, $self->length, 'color' => $par{'symcolor'};

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

	$c1 = $self->{'string'}->raw($i); $c2 = $othr->{'string'}->raw($i);

	#warn "[$i]= $c1 <=> $c2\n";

	#white space: no color
	next    if $self->{'string'}->is_space($c1);
					 
	#gap: gapcolour
	if ($self->{'string'}->is_non_sequence($c1)) {
	    push @$color, $i, 'color' => $par{'gapcolor'};
	    next;
	}

	#same symbol when upcased: use symbol/wildcard color
	if (uc $c1 eq uc $c2) {

            @tmp = $self->get_color($c1, $par{'colormap'});

	    if (@tmp) {
		if ($par{'css1'}) {
		    push @$color, $i, 'class' => $tmp[1];
		} else {
		    push @$color, $i, 'color' => $tmp[0];
		}
	    } else {
		push @$color, $i, 'color' => $par{'symcolor'};
	    }

	    next;
	}

	#different symbol: use contrast colour
	#push @$color, $i, 'color' => $par{'symcolor'};
	
	#different symbol: use wildcard colour
	@tmp = $self->get_color('*', $par{'colormap'});
	if (@tmp) {
	    if ($par{'css1'}) {
		push @$color, $i, 'class' => $tmp[1];
	    } else {
		push @$color, $i, 'color' => $tmp[0];
	    }
	} else {
	    push @$color, $i, 'color' => $par{'symcolor'};
	}
    }

    $self->{'display'}->{'paint'}  = 1;
    $self;
}

sub set_identity {
    #warn "Bio::MView::Align::Sequence::set_identity(@_)\n";
    my $self = shift;
    my $ref = shift;
    my $mode = shift;
    my $val = $self->compute_identity_to($ref, $mode);
    $self->set_display('label4'=>sprintf("%.1f%%", $val));
}

#compute percent identity to input reference object. Normalisation
#depends on the mode argument: 'reference' divides by the reference
#sequence length,'aligned' by the aligned region length (like blast),
#and 'hit' by the hit sequence. The last is the same as 'aligned' for
#blast, but different for multiple alignments like clustal.
sub compute_identity_to {
    #warn "Bio::MView::Align::Sequence::compute_identity_to(@_)\n";
    my $self = shift;
    my $othr = shift;
    my $mode = shift;
    my ($end, $i, $c1, $c2, $sum, $len, $gap);

    return unless $othr;

    die "${self}::compute_identity_to() length mismatch\n"
	unless $self->length == $othr->length;
    
    $end = $self->length +1;
    $sum = $len = 0;

    for ($i=1, $gap=0; $i<$end; $i++, $gap=0) {

	$c1 = $self->{'string'}->raw($i); $c2 = $othr->{'string'}->raw($i);
	
	$gap++  if $self->{'string'}->is_non_sequence($c1);
	$gap++  if $self->{'string'}->is_non_sequence($c2);
	
	next    if $gap > 1;
	
	#zero or at most one gap:

	#ignore leader/trailer gaps
	$len++  unless $self->{'string'}->is_padding($c1);

	$sum++  if $c1 eq $c2;
    }
    
    #normalise identities
    my $norm = 0;
    if ($mode eq 'reference') {
	$norm = $othr->ungapped_length;
	#warn "ref norm= $norm\n";
    } elsif ($mode eq 'aligned') {
	$norm = $len;
	#warn "aln norm= $norm\n";
    } elsif ($mode eq 'hit') {
	$norm = $self->ungapped_length;
	#warn "hit norm= $norm\n";
    }
    #warn "identity $self->{'id'} = $sum/$norm\n";

    return ($sum = 100 * ($sum + 0.0) / $norm)    if $norm > 0;
    return 0;
}

#compute percent identity to input reference object over length of aligned
#region of reference (same as blast). The latter is actually calculated by
#looking at the length of self, minus any terminal gaps.
sub find_identical_to {
    my $self = shift;
    my $othr = shift;
    my ($end, $i, $c1, $c2, $sum, $len, $gap, $s, $new);

    return unless $othr;

    die "${self}::find_identical_to() length mismatch\n"
	unless $self->length == $othr->length;
    
    $end = $self->length +1;
    $sum = $len = 0;
    $s   = '';

    for ($i=1, $gap=0; $i<$end; $i++, $gap=0) {

	$c1 = $self->{'string'}->raw($i); $c2 = $othr->{'string'}->raw($i);

	$gap++  if $self->{'string'}->is_non_sequence($c1);
	$gap++  if $self->{'string'}->is_non_sequence($c1);

	if ($gap > 1) {
	    $s .= $Bio::MView::Sequence::Text_Spc;
	    next;
	}

	#zero or at most one gap:

	#ignore leader/trailer gaps
	$len++  unless $self->{'string'}->is_padding($c1);

	if ($c1 ne $c2) {
	    $s .= $Bio::MView::Sequence::Text_Spc;
	    next;
	}

	$s .= $c1;
	$sum++;
    }

    $sum = 100 * ($sum + 0.0) / $len    if $len > 0;
   
    #encode the new "sequence"
    $new = new Bio::MView::Sequence;
    $new->append([\$s, $self->{'from'}, $self->{'from'}+$end-2]);

    #return a new identity Row therefrom
    Bio::MView::Align::Identity->new($self->id, $othr->id, $new, $sum);
}


###########################################################################
1;
