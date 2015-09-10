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

# $Id: Align.pm,v 1.47 2013/12/01 18:56:16 npb Exp $

######################################################################
package Bio::MView::Align;

use Bio::MView::Sequence;
use Bio::MView::Display;
use Bio::MView::Align::Row;

use strict;

use vars qw($Colormaps $Palette $Map_Text);

$Colormaps   = {};          #static hash of colormaps
$Palette     = [{},[]];     #static color palette
$Map_Text    = '';          #used as special index

my $BLOCKSEPARATOR = ':';   #for block search patterns

my %Template = 
    (
     'length'     => 0,     #alignment width
     'id2index'   => undef, #hash of identifiers giving row numbers
     'index2row'  => undef, #ordered list by row number of aligned objects
     'parent'     => undef, #identifier of of parent sequence
     'cursor'     => -1,    #index2row iterator
     'ref_id'     => undef, #identifier of reference row
     'tally'      => undef, #column tallies for consensus
     'coloring'   => undef, #coloring mode
     'colormap'   => undef, #name of colormap
     'colormap2'  => undef, #name of second colormap
     'group'      => undef, #consensus group name
     'ignore'     => undef, #ignore self/non-self classes
     'con_gaps'   => undef, #ignore gaps when computing consensus
     'threshold'  => undef, #consensus threshold for colouring
     'bold'       => undef, #display alignment in bold
     'css1'       => undef, #use CSS1 style sheets
     'alncolor'   => undef, #colour of alignment background
     'symcolor'   => undef, #default colour of alignment text
     'gapcolor'   => undef, #colour of alignment gap
     'old'        => {},    #previous settings of the above
     'nopshash'   => undef, #hash of id's to ignore for computations/colouring
     'hidehash'   => undef, #hash of id's to ignore for display
     'find'       => undef, #pattern to match in sequence
    );

my %Known_Parameter = 
    (
     'ref_id'     => [ '(\S+(?:\s*)?)+', undef ],
     'coloring'   => [ '\S+',     'none' ], 
     'colormap'   => [ '\S+',     $Bio::MView::Align::Sequence::Default_Colormap ],
     'colormap2'  => [ '\S+',     $Bio::MView::Align::Consensus::Default_Colormap ],
     'bold'       => [ '[01]',    1 ],
     'css1'       => [ '[01]',    0 ],
     'alncolor'   => [ '\S+',     $Bio::MView::Align::Row::Colour_White ],
     'labcolor'   => [ '\S+',     $Bio::MView::Align::Row::Colour_Black ],
     'symcolor'   => [ '\S+',     $Bio::MView::Align::Row::Colour_Black ],
     'gapcolor'   => [ '\S+',     $Bio::MView::Align::Row::Colour_DarkGray ],
     'group'      => [ '\S+',     $Bio::MView::Align::Consensus::Default_Group ],
     'ignore'     => [ '\S+',     $Bio::MView::Align::Consensus::Default_Ignore ],
     'con_gaps'   => [ '[01]',    1 ],
     'threshold'  => [ [],        [80] ],
     'nopshash'   => [ {},        {} ],
     'hidehash'   => [ {},        {} ],
     'find'       => [ '\S*',     '' ],
    );

my %Known_Molecule_Type =
    (
     #name
     'aa'         => 1,    #protein
     'na'         => 1,    #DNA/RNA
    );

my %Known_Alignment_Color_Schemes =
    (
     #name
     'none'       => 1,
     'any'        => 1,
     'identity'   => 1,
     'consensus'  => 1,
     'group'      => 1,
    );

my %Known_Consensus_Color_Schemes =
    (
     #name
     'none'       => 1,
     'any'        => 1,
     'identity'   => 1,
    );

#static load the $Colormaps hash.
eval { load_colormaps() }; if ($@) {$::COMPILE_ERROR=1; warn $@}


sub new {
    my $type = shift;
    #warn "${type}::new() @_\n";
    if (@_ < 1) {
	die "${type}::new() missing arguments\n";
    }
    my ($obj, $parent) = (@_, undef);
    my $i;

    my %self = %Template;

    $self{'id2index'}  = {};
    $self{'index2row'} = [];

    for ($i=0; $i<@$obj; $i++) {

	if (defined $obj->[$i]) {

	    #warn "[$i] ",  $obj->[$i]->id, " ", $obj->[$i]->string, "\n";

	    $self{'id2index'}->{$obj->[$i]->id} = $i;
	    $self{'index2row'}->[$i] = $obj->[$i];

	    $self{'length'} = $obj->[$i]->length    if $self{'length'} < 1;
	    
	    if ($obj->[$i]->length != $self{'length'}) {
		die "${type}::new() incompatible alignment lengths, row $i, expect $self{'length'}, got @{[$obj->[$i]->length]}\n";
	    }

	}
    }

    if (defined $parent) {
	$self{'parent'} = $parent;
    } else {
	$self{'parent'} = $self{'index2row'}->[0];
    } 

    my $self = bless \%self, $type;
    $self->initialise_parameters;
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
    foreach my $i (@{$self->{'index2row'}}) {
	$i->print if defined $i;
    }
    warn "\n";
    $self;
}

sub initialise_parameters {
    my $self = shift;
    my ($p) = (@_, \%Known_Parameter);
    local $_;
    foreach (keys %$p) {
	#warn "initialise_parameters() $_\n";
	if (ref $p->{$_}->[0] eq 'ARRAY') {
	    $self->{$_} = [];
	    next;
	}
	if (ref $p->{$_}->[0] eq 'HASH') {
	    $self->{$_} = {};
	    next;
	}
	$self->{$_} = $p->{$_}->[1];
    }

    $self->{'nopshash'} = {};
    $self->{'hidehash'} = {};

    $self;
}

sub set_parameters {
    my $self = shift;
    my $p = ref $_[0] ? shift : \%Known_Parameter;
    my ($key, $val);
    #warn "set_parameters($self) ". join(" ", keys %$p), "\n";
    while ($key = shift) {
	$val = shift;
	#warn "set_parameters() $key, $val\n";
	if (exists $p->{$key}) {
	    #warn "set_parameters() $key, $val\n";
	    if (ref $p->{$key}->[0] eq 'ARRAY' and ref $val eq 'ARRAY') {
		$self->{'old'}->{$key} = $self->{$key};
		$self->{$key} = $val;
		next;
	    }
	    if (ref $p->{$key}->[0] eq 'HASH' and ref $val eq 'HASH') {
		$self->{'old'}->{$key} = $self->{$key};
		$self->{$key} = $val;
		next;
	    }
	    if (! defined $val) {
		#set default
		$self->{'old'}->{$key} = $self->{$key};
		$self->{$key} = $p->{$key}->[1];
		next;
	    }
	    if ($val =~ /^$p->{$key}->[0]$/) {
		#matches expected format
		$self->{'old'}->{$key} = $self->{$key};
		$self->{$key} = $val;
		next;
	    }
	    warn "${self}::set_parameters() bad value for '$key', got '$val', wanted '$p->{$key}->[0]'\n";
	}
	#ignore unrecognised parameters which may be recognised by subclass
	#set_parameters() methods.
	#warn "set_parameters(IGNORE) $key, $val\n";
    }

    $self;
}

#concatenate Align object Rows into a new Align object: can be used to copy
sub cat {
    my $self = $_[0];
    my $i;
    my @obj = ();
    foreach (@_) {
	for ($i=0; $i<@{$_->{'index2row'}}; $i++) {
	    next unless defined $_->{'index2row'}->[$i];
	    push @obj, $_->{'index2row'}->[$i];
	}
    }
    new Bio::MView::Align(\@obj, $self->{'parent'});
}

sub length { $_[0]->{'length'} }

#return list of identifiers
sub ids {
    my @id = ();
    foreach (@{$_[0]->{'index2row'}}) {
	push @id, $_->id    if defined $_;
    }
    @id;
}

#return list of visible identifiers
sub visible { 
    my @id = ();
    foreach (@{$_[0]->{'index2row'}}) {
	push @id, $_->id
	    if defined $_ and ! exists $_[0]->{'hidehash'}->{$_->id};
    }
    @id;
}

#return number of stored rows
sub rows   { scalar map { $_->id if defined $_ } @{$_[0]->{'index2row'}} };

#return row object indexed by identifier
sub item {
    my ($self, $id) = @_;
    return 0    unless defined $id;
    if (exists $self->{'id2index'}->{$id}) {
	return $self->{'index2row'}->[$self->{'id2index'}->{$id}];
    }
    0;
}

#delete row(s) by identifier
sub delete {
    my $self = shift;
    local $_;
    foreach (@_) {
	$self->{'index2row'}->[$self->{'id2index'}->{$_}] = undef;
	$self->{'id2index'}->{$_} = undef;
    }
    $self;
}

#initialise stream of row objects
sub reset { $_[0]->{'cursor'} = -1 }

#return next row object in stream, or return 0 and reinitialise
sub next {
    $_[0]->{'cursor'}++;
    if (defined $_[0]->{'index2row'}->[$_[0]->{'cursor'}]) {
	return $_[0]->{'index2row'}->[$_[0]->{'cursor'}];
    }
    $_[0]->{'cursor'} = -1;
    return 0;
}

#load colormap data from stream
sub load_colormaps {
    my ($stream, $override) = (@_, \*DATA, 1);    #allow replacement maps
    my ($state, $map, $c1, $c2, $seethru, $color, $de, $mapignore) = (0);
    local $_;
    while (<$stream>) {

	#comments, blank lines
	if (/^\s*\#/ or /^\s*$/) {
	    next  if $state != 1;
	    $de .= $_;
	    next;
	}

	#map [name]
	if (/^\s*\[\s*(\S+)\s*\]/) {
	    $map = uc $1;
	    if (exists $Colormaps->{$map} and !$override) {
		$mapignore = 1;  #just for duration of this map
	    } else {
		$mapignore = 0;
	    }
	    $state = 1;
	    $de = '';
	    next;
	}

	#palette data: colorname: RGBcode
	#colorname MUST begin with a letter
	if (/^\s*colou?r\s+([a-z][-a-z0-9_]*)\s*:\s*(\#[0123456789ABCDEF]{6})/i) {
	    if (! exists $Palette->[0]->{$1}) {
		#warn "palette |$1|$2|\n";
		push @{$Palette->[1]}, $2;
		$Palette->[0]->{$1} = $#{$Palette->[1]};  #name->index
		$Palette->[0]->{$#{$Palette->[1]}} = $1;  #index->name
	    }
	    next;
	}

	die "load_colormaps(): colormap name undefined\n"
	    unless defined $map;

	next  if $mapignore;    #forget it if we're not allowing overrides

	#save map description?
	$Colormaps->{$map}->{$Map_Text} = $de  if $state == 1;
	
	#symbol[symbol] {->|=>} [palette_colorname|#RGB] any other text...
	if (/^\s*(\S)(\S)?\s*(->|=>)\s*(\S+)(?:\s+(.*))?/i) {

	    ($c1, $c2, $seethru, $color, $de) =
		($1, $2, ($3 eq '->' ? 'T' : 'S'), $4, (defined $5 ? $5 : ''));

	    $state = 2;

	    #only allow new colors in form of RGB codes
	    if (! exists $Palette->[0]->{$color}) {
		if ($color =~ /\#[0123456789ABCDEF]{6}/) {
		    #warn "new color  |$color| (transparency=$seethru)\n";
		    push @{$Palette->[1]}, $color;
		    $Palette->[0]->{$color} = $#{$Palette->[1]};  #name->index
		    $Palette->[0]->{$#{$Palette->[1]}} = $color;  #index->name
		} else {
		    chomp;
		    die "load_colormaps(): undefined color in line '$_'\n";
		}
	    }

	    $Colormaps->{$map}->{$c1} = [$Palette->[0]->{$color},$seethru,$de];
	    $Colormaps->{$map}->{$c2} = [$Palette->[0]->{$color},$seethru,$de]
		if defined $c2;

	    next;
	}

	#default
	chomp;
	die "Bio::MView::Align::load_colormaps(): bad format in line '$_'\n";
    }
    close $stream;
}

#return a descriptive listing of all known colormaps
sub list_colormaps {
    my $html = shift;
    my ($map, $sym, %p, $u, $l, $col, $rgb);
    my ($s, $f1, $f2, $c0, $c1, $c2) = ('', '', '', '', '', '');

    ($c0, $c1, $c2) = (
	"<FONT COLOR=\"$Bio::MView::Align::Row::Colour_Black\">",
	"<FONT COLOR=\"$Bio::MView::Align::Row::Colour_Comment\">",
	"</FONT>")  if $html;

    $s .= "$c1#Palette:\n\n#color                     : #RGB$c2\n";

    for ($u=0; $u < @{$Palette->[1]}; $u++) {
	$rgb = $Palette->[1]->[ $u ];
	($f1, $f2) = ("<FONT COLOR=\"$rgb\">", "</FONT>")  if $html;
	$s .= sprintf("color %s%-20s%s : %-7s\n", 
		      $f1, $Palette->[0]->{$u}, $f2, $Palette->[1]->[$u]);
    }

    $s .= "\n\n$c1#Colormap listing - suitable for reloading.\n";
    $s .= "#Character matching is case-sensitive.$c2\n\n";
    
    @_ = keys %$Colormaps  unless @_;

    foreach $map (sort @_) {
	$s .= "$c0\[$map]$c2\n";
	$s .= "$c1$Colormaps->{$map}->{$Map_Text}";
	$s .= "#symbols =>  color                #comment$c2\n";

	%p = %{$Colormaps->{$map}};    #copy colormap structure

	foreach $sym (sort keys %p) {

	    next    if $sym eq $Map_Text;
	    next    unless exists $p{$sym} and defined $p{$sym};

	    ($u, $l) = (uc $sym, lc $sym);

	    $col = $Palette->[0]->{ $p{$sym}->[0] };
	    $rgb = $Palette->[1]->[ $p{$sym}->[0] ];

	    ($f1, $f2) = ("<FONT COLOR=\"$rgb\">", "</FONT>")  if $html;

	    #lower and upper case: same symbol
	    if ($u eq $l) {
		$s .= sprintf("%s%-7s%s  %s  %s%-20s%s $c1%s$c2\n",
			      $f1, $sym, $f2,
			      ($p{$sym}->[1] eq 'T' ? '->' : '=>'),
			      $f1, $col, $f2,
			      $p{$sym}->[2]);
		next;
	    }

	    #lower and upper case: two symbols
	    if (exists $p{$u} and exists $p{$l}) {

                if ($p{$u}->[0] eq $p{$l}->[0] and
                    $p{$u}->[1] eq $p{$l}->[1]) {

		    #common definition
		    $s .= sprintf("%s%-7s%s  %s  %s%-20s%s $c1%s$c2\n",
				  $f1, $u . $l, $f2,
				  ($p{$sym}->[1] eq 'T' ? '->' : '=>'),
				  $f1, $col, $f2,
				  $p{$sym}->[2]);

		    $p{$u} = $p{$l} = undef;    #forget both

		} else {
		    #different definitions
		    $s .= sprintf("%s%-7s%s  %s  %s%-20s%s $c1%s$c2\n",
				  $f1, $sym, $f2,
				  ($p{$sym}->[1] eq 'T' ? '->' : '=>'),
				  $f1, $col, $f2,
				  $p{$sym}->[2]);
		}
		next;
	    }

	    #default: single symbol
	    $s .= sprintf("%s%-7s%s  %s  %s%-20s%s $c1%s$c2\n",
			  $f1, $sym, $f2,
			  ($p{$sym}->[1] eq 'T' ? '->' : '=>'),
			  $f1, $col, $f2,
			  $p{$sym}->[2]);
	    next;
	}

	$s .= "\n";
    }
    $s;
}

sub get_colormap_length {
    my $map = uc shift;
    die "get_colormap_length(): colormap name unknown\n"
	unless exists $Colormaps->{$map};
    my $len = scalar keys %{$Colormaps->{$map}};
    $len--  if defined $Colormaps->{$map}->{''};
    if ($len % 2 != 0) {
	die "get_colormap_length(): divide by two error\n"
    }
    return $len / 2;
}

sub print_css1_colormaps {
    my %color = @_;
    my ($s, $i, $rgb, $fg);

    $color{'alncolor'} = $Known_Parameter{'alncolor'}->[1]
        unless exists $color{'alncolor'};
    $color{'labcolor'} = $Known_Parameter{'labcolor'}->[1]
        unless exists $color{'labcolor'};
    $color{'symcolor'} = $Known_Parameter{'symcolor'}->[1]
        unless exists $color{'symcolor'};
    $color{'gapcolor'} = $Known_Parameter{'gapcolor'}->[1]
        unless exists $color{'gapcolor'};
    
    #warn "bg=$color{'alncolor'} fg=$color{'symcolor'}\n";
    
    $s = "TD{font-family:Fixed,Courier,monospace;background-color:$color{'alncolor'};color:$color{'labcolor'};}\n";

    for ($i=0; $i < @{$Palette->[1]}; $i++) {
	
	$rgb = $Palette->[1]->[$i];
	$fg = hex(substr($rgb, 1));
	my ($r, $g, $b) = (($fg>>16)&255, ($fg>>8)&255, $fg&255);
	
#	#SOLID: coloured background/monochrome foreground
#	#flip foreground between black/white depending on
#	#largest (brightest) RGB component of background.
#	$fg = $r; $fg = $g  if $g > $fg; $fg = $b  if $b > $fg;
#	if ($fg > 200) {
#	    $fg = $Bio::MView::Align::Row::Colour_Black;
#	} else {
#	    $fg = $Bio::MView::Align::Row::Colour_White;
#	}

	#just look at the green component
	if ($g > 200) {
	    $fg = $Bio::MView::Align::Row::Colour_Black;	    
	} else {
	    $fg = $Bio::MView::Align::Row::Colour_White;
	}

	#block: backgronud + foreground
	$s .= "FONT.S${i}\{background-color:$rgb;color:$fg} ";

	#transparent: no background + coloured foreground
	$s .= "FONT.T${i}\{color:$rgb} ";

	#comment
	$s .= "/* $Palette->[0]->{$i} */\n";
    }

    $s;
}

sub check_molecule_type {
    if (defined $_[0]) {
	if (exists $Known_Molecule_Type{lc $_[0]}) {
	    return lc $_[0];
	}
	return undef;
    }
    return map { lc $_ } sort keys %Known_Molecule_Type;
}

sub check_alignment_color_scheme {
    if (defined $_[0]) {
	if (exists $Known_Alignment_Color_Schemes{lc $_[0]}) {
	    return lc $_[0];
	}
	return undef;
    }
    return map { lc $_ } sort keys %Known_Alignment_Color_Schemes;
}

sub check_consensus_color_scheme {
    if (defined $_[0]) {
	if (exists $Known_Consensus_Color_Schemes{lc $_[0]}) {
	    return lc $_[0];
	}
	return undef;
    }
    return map { lc $_ } sort keys %Known_Consensus_Color_Schemes;
}

sub check_colormap {
    if (defined $_[0]) {
	if (exists $Colormaps->{uc $_[0]}) {
	    return uc $_[0];
	}
        if (exists $Palette->[0]->{lc $_[0]}) {
	    return lc $_[0]; #colormap all one predefined colour
        }
	return undef;
    }
    return sort keys %$Colormaps;
}

sub get_default_colormaps {
    if (! defined $_[0] or $_[0] eq 'aa') {
	#default to protein
	return ($Bio::MView::Align::Sequence::Default_PRO_Colormap,
		$Bio::MView::Align::Consensus::Default_PRO_Colormap);
    }
    #otherwise DNA/RNA explicitly requested
    return ($Bio::MView::Align::Sequence::Default_DNA_Colormap,
	    $Bio::MView::Align::Consensus::Default_DNA_Colormap);
}

sub get_default_find_colormap {
    return $Bio::MView::Align::Sequence::Default_FIND_Colormap;
}

#propagate display parameters to row objects
sub set_display {
    my $self = shift;
    local $_;
    foreach (@{$self->{'index2row'}}) {
	if (defined $_) {
	    $_->set_display(@_);
	}
    }
}

#ignore id's in remaining arglist
sub set_identity {
    #warn "Bio::MView::Align::set_identity(@_)\n";
    my $self = shift;
    my $ref  = shift;
    my $mode = shift;
    my $nops = shift;

    return unless defined $self->{'id2index'}->{$ref};
    return unless defined $self->{'index2row'}->[$self->{'id2index'}->{$ref}];

    $ref = $self->{'index2row'}->[$self->{'id2index'}->{$ref}];

    foreach (@{$self->{'index2row'}}) {
	if (defined $_ and ! exists $nops->{$_->id}) {
	    $_->set_identity($ref, $mode);
	}
    }
}

sub header {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s    if $quiet;

    if ($self->{'coloring'} eq 'any') {
	$s .= "Colored by: property\n";
    }
    elsif ($self->{'coloring'} eq 'identity' and defined $self->{'ref_id'}) {
	$s .= "Colored by: identity + property\n";
    }
    elsif ($self->{'coloring'} eq 'consensus') {
	$s .= "Colored by: consensus/$self->{'threshold'}->[0]\% and property\n";
    }
    elsif ($self->{'coloring'} eq 'group') {
	$s .= "Colored by: consensus/$self->{'threshold'}->[0]\% and group property\n";
    }
    elsif ($self->{'coloring'} eq 'find') {
	$s .= "Colored by: search pattern '$self->{'find'}'\n";
    }
    Bio::MView::Display::displaytext($s);
}

sub set_color_scheme {
    my $self = shift;

    $self->set_parameters(@_);

    return $self    if $self->{'coloring'} eq 'none';

    #user-defined colouring?

    $self->color_special('colormap'  => $self->{'colormap'},
			 'colormap2' => $self->{'colormap2'},
			 'symcolor'  => $self->{'symcolor'},
			 'gapcolor'  => $self->{'gapcolor'},
			 'css1'      => $self->{'css1'},
			);

    if ($self->{'coloring'} eq 'any') {
	$self->color_by_type('colormap'  => $self->{'colormap'},
			     'colormap2' => $self->{'colormap2'},
			     'symcolor'  => $self->{'symcolor'},
			     'gapcolor'  => $self->{'gapcolor'},
			     'css1'      => $self->{'css1'},
			    );
	return $self;
    }

    if ($self->{'coloring'} eq 'identity') {
	$self->color_by_identity($self->{'ref_id'},
				 'colormap'  => $self->{'colormap'},
				 'colormap2' => $self->{'colormap2'},
				 'symcolor'  => $self->{'symcolor'},
				 'gapcolor'  => $self->{'gapcolor'},
				 'css1'      => $self->{'css1'},
				);
	return $self;
    }

    if ($self->{'coloring'} eq 'consensus') {
	$self->color_by_consensus_sequence('colormap'  => $self->{'colormap'},
					   'colormap2' => $self->{'colormap2'},
					   'group'     => $self->{'group'},
					   'threshold' => $self->{'threshold'},
					   'symcolor'  => $self->{'symcolor'},
					   'gapcolor'  => $self->{'gapcolor'},
					   'css1'      => $self->{'css1'},
					  );
	return $self;
    }

    if ($self->{'coloring'} eq 'group') {
	$self->color_by_consensus_group('colormap'  => $self->{'colormap'},
					'colormap2' => $self->{'colormap2'},
					'group'     => $self->{'group'},
					'threshold' => $self->{'threshold'},
					'symcolor'  => $self->{'symcolor'},
					'gapcolor'  => $self->{'gapcolor'},
					'css1'      => $self->{'css1'},
				       );
	return $self;
    }

    if ($self->{'coloring'} eq 'find') {
	$self->color_by_find_block('colormap'  => $self->{'colormap'},
				   'colormap2' => $self->{'colormap2'},
				   'symcolor'  => $self->{'symcolor'},
				   'gapcolor'  => $self->{'gapcolor'},
				   'css1'      => $self->{'css1'},
				   'find'      => $self->{'find'},
			    );
	return $self;
    }

    warn "${self}::set_color_scheme() unknown mode '$self->{'coloring'}'\n";

    $self;
}

#propagate colour scheme to row objects
sub color_special {
    my $self = shift;

    for (my $i=0; $i<@{$self->{'index2row'}}; $i++) {
	next unless defined $self->{'index2row'}->[$i];
	next if exists $self->{'nopshash'}->{$self->{'index2row'}->[$i]->id};
	next if exists $self->{'hidehash'}->{$self->{'index2row'}->[$i]->id};
	next if $self->{'index2row'}->[$i]->{'type'} ne 'special';
	$self->{'index2row'}->[$i]->color_special(@_);
	$self->{'index2row'}->[$i]->set_display('label0'=>'',
						'label2'=>'',
						'label3'=>'',
						'label4'=>'',
						'label5'=>'',
						'label6'=>'',
						);
    }
    $self;
}

#propagate colour scheme to row objects
sub color_by_find_block {
    my $self = shift;
    my %par = @_;

    my $mapsize = get_colormap_length($par{'colormap'});
    my @patterns = split($BLOCKSEPARATOR, $par{find});
    if (@patterns > $mapsize) {
	warn "recycling colormap '$par{colormap}': @{[scalar @patterns]} patterns but only $mapsize color(s)\n";
    }
    push @_, 'mapsize'  => $mapsize;
    push @_, 'patterns' => [@patterns];

    for (my $i=0; $i<@{$self->{'index2row'}}; $i++) {
	next unless defined $self->{'index2row'}->[$i];
	next if exists $self->{'nopshash'}->{$self->{'index2row'}->[$i]->id};
	next if exists $self->{'hidehash'}->{$self->{'index2row'}->[$i]->id};
	$self->{'index2row'}->[$i]->color_by_find_block(@_);
    }
    $self;
}

#propagate colour scheme to row objects
sub color_by_type {
    my $self = shift;
    
    for (my $i=0; $i<@{$self->{'index2row'}}; $i++) {
	next unless defined $self->{'index2row'}->[$i];
	next if exists $self->{'nopshash'}->{$self->{'index2row'}->[$i]->id};
	next if exists $self->{'hidehash'}->{$self->{'index2row'}->[$i]->id};
	$self->{'index2row'}->[$i]->color_by_type(@_);
    }
    $self;
}

#propagate colour scheme to row objects
sub color_by_identity {
    my ($self, $id) = (shift, shift);
    my $ref = $self->item($id);

    for (my $i=0; $i<@{$self->{'index2row'}}; $i++) {
	next unless defined $self->{'index2row'}->[$i];
	next if exists $self->{'nopshash'}->{$self->{'index2row'}->[$i]->id};
	next if exists $self->{'hidehash'}->{$self->{'index2row'}->[$i]->id};
	$self->{'index2row'}->[$i]->color_by_identity($ref, @_);
    }
    $self;
}

#propagate colour scheme to row objects
sub color_by_consensus_sequence {
    my $self = shift;

    #is there already a suitable tally?
    if (!defined $self->{'tally'} or
	(defined $self->{'old'}->{'group'} and 
	 $self->{'old'}->{'group'} ne $self->{'group'})) {
	$self->compute_tallies($self->{'group'});
    }

    my $con = new Bio::MView::Align::Consensus($self->{'tally'},
					       $self->{'group'},
					       $self->{'threshold'}->[0],
					       $self->{'ignore'},
					       $self->{'parent'}->{'from'});
    
    for (my $i=0; $i<@{$self->{'index2row'}}; $i++) {
	next unless defined $self->{'index2row'}->[$i];
	next if exists $self->{'nopshash'}->{$self->{'index2row'}->[$i]->id};
	next if exists $self->{'hidehash'}->{$self->{'index2row'}->[$i]->id};
	$con->color_by_consensus_sequence($self->{'index2row'}->[$i], @_);
    }
    $self;
}

#propagate colour scheme to row objects
sub color_by_consensus_group {
    my $self = shift;

    #is there already a suitable tally?
    if (!defined $self->{'tally'} or
	(defined $self->{'old'}->{'group'} and 
	 $self->{'old'}->{'group'} ne $self->{'group'})) {

	$self->compute_tallies($self->{'group'});
    }

    my $con = new Bio::MView::Align::Consensus($self->{'tally'},
					       $self->{'group'},
					       $self->{'threshold'}->[0],
					       $self->{'ignore'},
					       $self->{'parent'}->{'from'});
    
    for (my $i=0; $i<@{$self->{'index2row'}}; $i++) {
	next unless defined $self->{'index2row'}->[$i];
	next if exists $self->{'nopshash'}->{$self->{'index2row'}->[$i]->id};
	next if exists $self->{'hidehash'}->{$self->{'index2row'}->[$i]->id};
	$con->color_by_consensus_group($self->{'index2row'}->[$i], @_);
    }
    $self;
}

#return array of Bio::MView::Display::display() constructor arguments
sub init_display { ( $_[0]->{'parent'}->{'string'} ) }

#append Row data to the input Display object: done one at a time to 
#reduce memory usage instead of accumulating a potentially long list before
#passing to Display::append(), and to permit incremental garbage collection of
#each Align::Row object oce it has been appended. the latter can be switched
#off if optional argument $nogc is true (needed when further processing of Row
#objects will occur, eg., consensus calculations).
sub append_display {
    my ($self, $dis, $nogc) = (@_, 0);
    #warn "append_display($dis, $nogc)\n";
    for (my $i=0; $i<@{$self->{'index2row'}}; $i++) {
	if (defined $self->{'index2row'}->[$i]) {
	    next  if
		exists $self->{'hidehash'}->{$self->{'index2row'}->[$i]->id};
	    
	    #append the row data structure to the Display object
	    $dis->append($self->{'index2row'}->[$i]->get_display);

	    #optional garbage collection
	    $self->{'index2row'}->[$i] = undef  unless $nogc;
	}
    }
    $self;
}

#compute effective all pairwise alignment and keep only those sequences
#with $min <= pairwise identity <= $max. also keep any sequences with id's
#supplied as remaining arguments. The mode argument determines which
#sequence length is used to normalise: 'reference', 'aligned', 'hit'.
sub prune_all_identities {
    my ($self, $mode, $min, $max, $topn) = (shift, shift, shift, shift, shift);
    my (@obj, %keep, $ref, $i, $row);

    @obj = ();
    $min = 0    if $min < 0;
    $max = 100  if $max > 100;

    #special case
    return $self  unless $min > 0 or $max < 100;
    #return $self  if $min > $max;  #silly combination

    #ensure no replicates in keep list
    foreach $i (@_) {
	$ref = $self->{'index2row'}->[$self->{'id2index'}->{$i}];
	$keep{$ref} = $ref    if defined $ref;
    }

    #prime keep list
    @obj = ();

    #compare all rows not on keep list against latter and add thereto if
    #sufficiently dissimilar
    for ($i=0; $i<@{$self->{'index2row'}}; $i++) {
	
	next unless defined $self->{'index2row'}->[$i];

	#enforce limit on number of rows
	last    if $topn > 0 and @obj == $topn;

	if (exists $keep{$self->{'index2row'}->[$i]}) {
	    push @obj, $self->{'index2row'}->[$i];
	    next;
	}

	$row = $self->{'index2row'}->[$i];
	
	foreach $ref (@obj) {
	    my $pcid = $row->compute_identity_to($ref, $mode);
	    #store object if %identity satisfies cutoff for all kept hits
	    if ($pcid < $min or $pcid > $max) {
		$row = 0;
		last;
	    }
	}
	
	if ($row) {
	    #print STDERR "passed ", $row->id, "\n";
	    push @obj, $row;
	}

	#warn join(" ", map { $_->id } @obj), "\n";
    }

    new Bio::MView::Align(\@obj, $self->{'parent'});
}

#generate a new alignment from an existing one with extra information
#showing %identities and identical symbols with respect to some supplied
#identifier. only keep lines  with %identity to reference <= $limit.
sub prune_identities_gt {
    my ($self, $id, $limit, $mode) = @_;
    my ($ref, $i, @obj, $row, $val);
    
    $ref = $self->item($id);

    @obj = ();

    for ($i=0; $i<@{$self->{'index2row'}}; $i++) {

	next unless defined $self->{'index2row'}->[$i];

	$row = $self->{'index2row'}->[$i];

	#store object if %identity satisfies cutoff OR if the object was
	#the reference object!
	if (($val = $row->compute_identity_to($ref, $mode)) <= $limit or
	     $row->id eq $id) {

	    $row->set_display('label4'=>sprintf("%.1f%%", $val));
	
	    push @obj, $row;
	}
    }

    new Bio::MView::Align(\@obj, $self->{'parent'});
}

#generate a new alignment comprising a ruler based on this alignment
sub build_ruler {
    new Bio::MView::Align([new Bio::MView::Align::Ruler($_[0]->length)],
			  $_[0]->{'parent'});
}

#generate a new alignment using an existing one but with lines showing
#consensus sequences at specified percent thresholds
sub build_consensus_rows {
    my ($self, $group, $threshold, $ignore, $con_gaps) = @_;
    my ($thresh, $con, $i);
    
    $self->set_parameters('group' => $group, 'ignore' => $ignore,
			  'con_gaps' => $con_gaps);

    #is there already a suitable tally?
    if (!defined $self->{'tally'} or
	(defined $self->{'old'}->{'group'} and 
	 $self->{'old'}->{'group'} ne $self->{'group'})) {

	$self->compute_tallies($self->{'group'});
    }

    my @obj = ();
    
    foreach $thresh (@$threshold) {

	$con = new Bio::MView::Align::Consensus($self->{'tally'},
						$group, $thresh, $ignore,
						$self->{'parent'}->{'from'});

	$con->set_display('label0'=>'',
			  'label2'=>'',
			  'label3'=>'',
			  'label4'=>'',
			  'label5'=>'',
			  'label6'=>'',
			  );
	
	push @obj, $con;
    }

    new Bio::MView::Align(\@obj, $self->{'parent'});
}

sub compute_tallies {
    my ($self, $group) = @_;
    my ($row, $col, $r, $c);

    $group = $Bio::MView::Align::Consensus::Default_Group
	unless defined $group;

    $self->{'tally'} = [];

    #iterate over columns
    for ($c=1; $c <= $self->{'length'}; $c++) {

	$col = [];

	#iterate over rows
	for ($r=0; $r<@{$self->{'index2row'}}; $r++) {

	    $row = $self->{'index2row'}->[$r];
	    
	    next unless defined $row;
	    next if exists $self->{'nopshash'}->{$row->id};
	    next if exists $self->{'hidehash'}->{$row->id};
	    next if $row->{'type'} ne 'sequence';

	    push @$col, $row->{'string'}->raw($c);
	}

	#warn "compute_tallies: @$col\n";

	push @{$self->{'tally'}},
	    Bio::MView::Align::Consensus::tally($group, $col,
						$self->{'con_gaps'});
    }
    $self;
}


######################################################################
1;

__DATA__

######################################################################
# First some standard colours
######################################################################

# Netscape 216 (supposedly) cross-platform colours
# (other colours may not resolve well on all platforms)

#Primaries/secondaries (cross-platform colours)
color  black               :  #000000
color  white               :  #ffffff
color  red                 :  #ff0000
color  green               :  #00ff00
color  blue                :  #0000ff
color  cyan                :  #00ffff
color  magenta             :  #ff00ff
color  yellow              :  #ffff00

#Miscellaneous (cross-platform colours)
color  purple              :  #6600cc
color  dull-blue           :  #0099ff
color  dark-green-blue     :  #33cccc
color  medium-green-blue   :  #00ffcc
color  bright-blue         :  #0033ff
color  dark-green          :  #009900
color  bright-green        :  #33cc00
color  orange              :  #ff3333
color  orange-brown        :  #cc6600
color  bright-red          :  #cc0000
color  light-gray          :  #999999
color  dark-gray           :  #666666

#Grey levels (some cross-platform)
color  gray0               :  #ffffff
color  gray1               :  #eeeeee
color  gray2               :  #dddddd
color  gray3               :  #cccccc
color  gray4               :  #bbbbbb
color  gray5               :  #aaaaaa
color  gray6               :  #999999
color  gray7               :  #888888
color  gray8               :  #777777
color  gray9               :  #666666
color  gray10              :  #555555
color  gray11              :  #444444
color  gray12              :  #333333
color  gray13              :  #222222
color  gray14              :  #111111
color  gray15              :  #000000

#CLUSTALX screen colours (protein + nucleotide)
color  clustal-red         :  #e53319
color  clustal-blue        :  #197fe5
color  clustal-green       :  #19cc19
color  clustal-cyan        :  #19b2b2
color  clustal-pink        :  #e57f7f
color  clustal-magenta     :  #cc4ccc
color  clustal-yellow      :  #cccc00
color  clustal-orange      :  #e5994c
color  clustal-white       :  #ffffff
color  clustal-light-gray  :  #999999 
color  clustal-dark-gray   :  #666666
color  clustal-black       :  #000000

#CLUSTALX printing colours (protein + nucleotide)
color  clustal-white-print 	:  #ffffff
color  clustal-yellow-print	:  #ffff00
color  clustal-violet-print     :  #6619e5
color  clustal-red-print   	:  #e57f66
color  clustal-blue-print  	:  #66e5e5
color  clustal-purple-print	:  #b299e5
color  clustal-black-print      :  #000000
color  clustal-pink-print       :  #e57f7f #overrides #cc4ccc in colprint.xml
color  clustal-cyan-print       :  #19b2b2
color  clustal-magenta-print    :  #cc4ccc
color  clustal-orange-print     :  #e5994c #overrides #e5b24c in colprint.xml
color  clustal-green-print      :  #19cc19
color  clustal-light-gray-print :  #999999
color  clustal-dark-gray-print  :  #99b2b2

#Frederico Nardi's suggested colours for CLUSTAL
color  nardi-red           :  #ff1111
color  nardi-blue          :  #1155ff
color  nardi-green         :  #11dd11
color  nardi-cyan          :  #11ffff
color  nardi-yellow        :  #ffff11
color  nardi-orange        :  #ff7f11
color  nardi-pink          :  #ff11ff
color  nardi-purple        :  #6611cc
color  nardi-dull-blue     :  #197fe5
color  nardi-dark-gray     :  #666666
color  nardi-light-gray    :  #999999 

#Kuang Lin's colour neural net derived scheme
color  lin-A               :  #90fe23
color  lin-R               :  #fe5e2d
color  lin-N               :  #2e3d2d
color  lin-D               :  #00903b
color  lin-C               :  #004baa
color  lin-Q               :  #864b00
color  lin-E               :  #3fa201
color  lin-G               :  #10fe68
color  lin-H               :  #b2063b
color  lin-I               :  #04ced9
color  lin-L               :  #4972fe
color  lin-K               :  #c4a100
color  lin-M               :  #2a84dd
color  lin-F               :  #a60ade
color  lin-P               :  #fe61fe
color  lin-S               :  #f7e847
color  lin-T               :  #fefeb3
color  lin-W               :  #4a007f
color  lin-Y               :  #e903a8
color  lin-V               :  #5bfdfd

#block colours
color  find-A             :  #90fe23
color  find-B             :  #fe5e2d
color  find-C             :  #2e3d2d
color  find-D             :  #00903b
color  find-E             :  #004baa
color  find-F             :  #864b00
color  find-G             :  #3fa201
color  find-H             :  #10fe68
color  find-I             :  #b2063b
color  find-J             :  #04ced9
color  find-K             :  #4972fe
color  find-L             :  #c4a100
color  find-M             :  #2a84dd
color  find-N             :  #a60ade
color  find-O             :  #fe61fe
color  find-P             :  #f7e847
color  find-Q             :  #fefeb3
color  find-R             :  #4a007f
color  find-S             :  #e903a8
color  find-T             :  #5bfdfd


######################################################################
#now some colour schemes
######################################################################

#symbol -> colour (RGB hex or colorname)  [#comment]

[P1]
#Protein: highlight amino acid physicochemical properties
Gg  =>  bright-green         #hydrophobic
Aa  =>  bright-green         #hydrophobic
Ii  =>  bright-green         #hydrophobic
Vv  =>  bright-green         #hydrophobic
Ll  =>  bright-green         #hydrophobic
Mm  =>  bright-green         #hydrophobic
Ff  =>  dark-green           #large hydrophobic
Yy  =>  dark-green           #large hydrophobic
Ww  =>  dark-green           #large hydrophobic
Hh  =>  dark-green           #large hydrophobic
Cc  =>  yellow               #cysteine
Pp  =>  bright-green         #hydrophobic
Kk  =>  bright-red           #positive charge
Rr  =>  bright-red           #positive charge
Dd  =>  bright-blue          #negative charge
Ee  =>  bright-blue          #negative charge
Qq  =>  purple               #polar
Nn  =>  purple               #polar
Ss  =>  dull-blue            #small alcohol
Tt  =>  dull-blue            #small alcohol
Bb  =>  dark-gray            #D or N
Zz  =>  dark-gray            #E or Q
Xx  ->  dark-gray            #any
?   ->  light-gray           #unknown
*   ->  dark-gray            #mismatch

[GPCR]
#Protein: GPCRdb color scheme for Gert Vriend
Gg  =>  orange-brown         #backbone change
Pp  =>  orange-brown         #backbone change
Aa  =>  bright-green         #hydrophobic
Ii  =>  bright-green         #hydrophobic
Vv  =>  bright-green         #hydrophobic
Ll  =>  bright-green         #hydrophobic
Mm  =>  bright-green         #hydrophobic
Cc  =>  yellow               #cysteine
Qq  =>  bright-red           #positive charge
Ee  =>  bright-red           #positive charge
Nn  =>  bright-red           #positive charge
Dd  =>  bright-red           #positive charge
Bb  =>  bright-red           #D or N
Zz  =>  bright-red           #E or Q
Hh  =>  bright-blue          #negative charge
Kk  =>  bright-blue          #negative charge
Rr  =>  bright-blue          #negative charge
Ss  =>  dark-green-blue      #small alcohol
Tt  =>  dark-green-blue      #small alcohol
Yy  =>  medium-green-blue    #large hydrophobic
Ff  =>  cyan                 #large hydrophobic
Ww  =>  cyan                 #large hydrophobic
Xx  ->  dark-gray            #any
?   ->  light-gray           #unknown
*   ->  dark-gray            #mismatch

[CYS]
#Protein: highlight cysteines
Cc  =>  yellow               #cysteine
Xx  ->  dark-gray            #any
?   ->  light-gray           #unknown
*   ->  dark-gray            #mismatch

[CHARGE]
#Protein: highlight charged amino acids
Kk  =>  bright-red           #positive charge
Rr  =>  bright-red           #positive charge
Dd  =>  bright-blue          #negative charge
Ee  =>  bright-blue          #negative charge
Bb  =>  dark-gray            #D or N
Zz  =>  dark-gray            #E or Q
Xx  =>  dark-gray            #any
?   =>  light-gray           #unknown
*   =>  black                #mismatch

[POLAR1]
#Protein: highlight charged and polar amino acids
Kk  =>  bright-red           #positive charge
Rr  =>  bright-red           #positive charge
Dd  =>  bright-blue          #negative charge
Ee  =>  bright-blue          #negative charge
Qq  =>  purple               #charged/polar
Nn  =>  purple               #charged/polar
Ss  =>  purple               #charged/polar
Tt  =>  purple               #charged/polar
Hh  =>  purple               #charged/polar
Bb  =>  purple               #D or N
Zz  =>  purple               #E or Q
Xx  ->  dark-gray            #any
?   ->  light-gray           #unknown
*   ->  dark-gray            #mismatch

[D1]
#DNA: highlight nucleotide types
Aa  =>  bright-blue          #purine
Gg  =>  bright-blue          #purine
Tt  =>  dull-blue            #pyrimidine
Cc  =>  dull-blue            #pyrimidine
Uu  =>  dull-blue            #pyrimidine
Mm  =>  dark-gray            #A or C
Rr  =>  dark-gray            #A or G
Ww  =>  dark-gray            #A or T
Ss  =>  dark-gray            #C or G
Yy  =>  dark-gray            #C or T
Kk  =>  dark-gray            #G or T
Vv  =>  dark-gray            #A or C or G; not T
Hh  =>  dark-gray            #A or C or T; not G
Dd  =>  dark-gray            #A or G or T; not C
Bb  =>  dark-gray            #C or G or T; not A
Nn  =>  dark-gray            #A or C or G or T
Xx  ->  dark-gray            #any
?   ->  light-gray           #unknown
*   ->  dark-gray            #mismatch

[D2]
#DNA: highlight match versus mismatch under consensus coloring schemes
Aa  =>  bright-blue          #match
Bb  =>  bright-blue          #match
Cc  =>  bright-blue          #match
Dd  =>  bright-blue          #match
Gg  =>  bright-blue          #match
Hh  =>  bright-blue          #match
Kk  =>  bright-blue          #match
Mm  =>  bright-blue          #match
Nn  =>  bright-blue          #match
Rr  =>  bright-blue          #match
Ss  =>  bright-blue          #match
Tt  =>  bright-blue          #match
Uu  =>  bright-blue          #match
Vv  =>  bright-blue          #match
Ww  =>  bright-blue          #match
Xx  =>  bright-blue          #match
Yy  =>  bright-blue          #match
?   ->  light-gray           #unknown
*   =>  bright-red           #mismatch

[CLUSTAL_NUC]
#CLUSTAL-derived colours for nucleotides
Aa  =>  clustal-red          #match
Bb  =>  clustal-dark-gray    #match
Cc  =>  clustal-blue         #match
Dd  =>  clustal-dark-gray    #match
Gg  =>  clustal-orange       #match
Hh  =>  clustal-dark-gray    #match
Kk  =>  clustal-dark-gray    #match
Mm  =>  clustal-dark-gray    #match
Nn  =>  clustal-dark-gray    #match
Rr  =>  clustal-dark-gray    #match
Ss  =>  clustal-dark-gray    #match
Tt  =>  clustal-green        #match
Uu  =>  clustal-green        #match
Vv  =>  clustal-dark-gray    #match
Ww  =>  clustal-dark-gray    #match
Xx  =>  clustal-dark-gray    #match
Yy  =>  clustal-dark-gray    #match
?   ->  clustal-light-gray   #unknown
*   =>  clustal-dark-gray    #mismatch

[PC1]
#Protein consensus: highlight equivalence class
a  ->  dark-green            #aromatic
l  ->  bright-green          #aliphatic
h  ->  bright-green          #hydrophobic
+  ->  bright-red            #positive charge
-  ->  bright-blue           #negative charge
c  ->  purple                #charged
p  ->  dull-blue             #polar
o  ->  dull-blue             #alcohol
u  ->  bright-green          #tiny
s  ->  bright-green          #small
t  ->  bright-green          #turnlike
*  ->  dark-gray             #mismatch

[DC1]
#DNA consensus: highlight ring type
r  ->  purple                #purine
y  ->  orange                #pyrimidine
*  ->  dark-gray             #mismatch

#+ Frederico Nardi
[NARDI]
#Protein: highlight amino acid physicochemical properties
Aa       =>  nardi-dull-blue    #hydrophobic
Bb       =>  nardi-pink         #D or N
Cc       =>  nardi-yellow       #cysteine
Dd       =>  nardi-pink         #negative charge
Ee       =>  nardi-pink         #negative charge
Ff       =>  nardi-dull-blue    #large hydrophobic
Gg       =>  nardi-orange       #glycine
Hh       =>  nardi-dull-blue    #large hydrophobic
Ii       =>  nardi-dull-blue    #hydrophobic
Kk       =>  nardi-red          #positive charge
Ll       =>  nardi-dull-blue    #hydrophobic
Mm       =>  nardi-dull-blue    #hydrophobic
Nn       =>  nardi-green        #polar
Pp       =>  nardi-yellow       #proline
Qq       =>  nardi-green        #polar
Rr       =>  nardi-red          #positive charge
Ss       =>  nardi-green        #small alcohol
Tt       =>  nardi-green        #small alcohol
Vv       =>  nardi-dull-blue    #hydrophobic
Ww       =>  nardi-dull-blue    #large hydrophobic
Xx       =>  nardi-dark-gray    #any
Yy       =>  nardi-cyan         #large hydrophobic
Zz       =>  nardi-pink         #E or Q
?        =>  nardi-light-gray   #unknown
*        =>  nardi-dark-gray    #mismatch

[CLUSTAL]
#Protein: highlight amino acid physicochemical properties
Aa       =>  clustal-blue         #hydrophobic
Bb       =>  clustal-white        #D or N
Cc       =>  clustal-pink         #hydrophobic
Dd       =>  clustal-magenta      #negative charge
Ee       =>  clustal-magenta      #negative charge
Ff       =>  clustal-blue         #large hydrophobic
Gg       =>  clustal-orange       #glycine
Hh       =>  clustal-cyan         #large hydrophobic
Ii       =>  clustal-blue         #hydrophobic
Kk       =>  clustal-red          #positive charge
Ll       =>  clustal-blue         #hydrophobic
Mm       =>  clustal-blue         #hydrophobic
Nn       =>  clustal-green        #polar
Pp       =>  clustal-yellow       #proline
Qq       =>  clustal-green        #polar
Rr       =>  clustal-red          #positive charge
Ss       =>  clustal-green        #small alcohol
Tt       =>  clustal-green        #small alcohol
Vv       =>  clustal-blue         #hydrophobic
Ww       =>  clustal-blue         #large hydrophobic
Xx       =>  clustal-dark-gray    #any
Yy       =>  clustal-cyan         #large hydrophobic
Zz       =>  clustal-white        #E or Q
?        =>  clustal-light-gray   #unknown
*        =>  clustal-dark-gray    #mismatch

[CCLUSTAL]
#Protein consensus: highlight equivalence class
+        ->  nardi-red          #positive charge
-        ->  nardi-pink         #negative charge
a        ->  nardi-dull-blue    #aromatic
c        ->  nardi-purple       #charged
h        ->  nardi-dull-blue    #hydrophobic
l        ->  nardi-dull-blue    #aliphatic
o        ->  nardi-green        #alcohol
p        ->  nardi-green        #polar
s        ->  nardi-blue         #small
t        ->  nardi-blue         #turnlike
u        ->  nardi-blue         #tiny
*        ->  nardi-dark-gray    #mismatch

[KXLIN]
#Kuang Lin's colour neural net derived scheme
Aa       =>  lin-A
Rr       =>  lin-R
Nn       =>  lin-N
Dd       =>  lin-D
Cc       =>  lin-C
Qq       =>  lin-Q
Ee       =>  lin-E
Gg       =>  lin-G
Hh       =>  lin-H
Ii       =>  lin-I
Ll       =>  lin-L
Kk       =>  lin-K
Mm       =>  lin-M
Ff       =>  lin-F
Pp       =>  lin-P
Ss       =>  lin-S
Tt       =>  lin-T
Ww       =>  lin-W
Yy       =>  lin-Y
Vv       =>  lin-V
Bb       =>  dark-gray            #D or N
Zz       =>  dark-gray            #E or Q
Xx       ->  dark-gray            #any
?        ->  light-gray           #unknown
*        ->  dark-gray            #mismatch

[FIND]
#colour successive 'find' pattern blocks
Aa       =>  find-A
Bb       =>  find-B
Cc       =>  find-C
Dd       =>  find-D
Ee       =>  find-E
Ff       =>  find-F
Gg       =>  find-G
Hh       =>  find-H
Ii       =>  find-I
Jj       =>  find-J
Kk       =>  find-K
Ll       =>  find-L
Mm       =>  find-M
Nn       =>  find-N
Oo       =>  find-O
Pp       =>  find-P
Qq       =>  find-Q
Rr       =>  find-R
Ss       =>  find-S
Tt       =>  find-T

##########################################################################
