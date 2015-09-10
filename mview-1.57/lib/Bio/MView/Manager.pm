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

# $Id: Manager.pm,v 1.55 2013/09/27 00:54:23 npb Exp $

######################################################################
package Bio::MView::Manager;

use Bio::MView::Build;
use Bio::MView::Display;
use Bio::SRS 'srsLink';
use strict;

my %Template =
    (
     'file'          => undef,
     'format'        => undef,
     'stream'        => undef,
     'filter'        => undef,
     'class'         => undef,
     'display'       => undef,
     'html'          => undef,
     'linkcolor'     => undef,
     'alinkcolor'    => undef,
     'vlinkcolor'    => undef,
     'register'      => undef,
     'quiet'         => undef,
     'noparse'       => undef,
     'bp'            => undef,
     'ap'            => undef,
    );

my %Template_Build_Param =
    (
     'topn'          => undef,
     'minident'      => undef,
     'maxident'      => undef,
     'pcid'          => undef,
     'mode'          => undef,
     'ref_id'        => undef,
     'keeplist'      => undef,
     'disclist'      => undef,
     'nopslist'      => undef,
     'range'         => undef,
     'gap'           => undef,
     		     
     'ruler'         => undef,
     'alignment'     => undef,
     'consensus'     => undef,
   		     
     'maxpval'       => undef,
     'maxeval'       => undef,
     'minbits'       => undef,
     'minscore'      => undef,
     'minopt'        => undef,
     'hsp'           => undef,
     		     
     'cycle'         => undef,
     'strand'        => undef,
     'chain'         => undef,
     'block'         => undef,

     'showpcid'      => undef,
    );

my %Template_Align_Param =
    (
     'aln_coloring'  => undef,
     'aln_colormap'  => undef,
     'aln_groupmap'  => undef,
     'aln_threshold' => undef,
     'aln_ignore'    => undef,
  		     
     'con_coloring'  => undef,
     'con_colormap'  => undef,
     'con_groupmap'  => undef,
     'con_threshold' => undef,
     'con_ignore'    => undef,
     'con_gaps'      => undef,
		     
     'css1'          => undef,
     'alncolor'      => undef,
     'labcolor'      => undef,
     'symcolor'      => undef,
     'gapcolor'      => undef,
     'bold'          => undef,
     'width'         => undef,
     'label0'        => undef,
     'label1'        => undef,
     'label2'        => undef,
     'label3'        => undef,
     'label4'        => undef,
     'label5'        => undef,
     'label6'        => undef,

     'find'          => undef,
    );

sub new {
    my $type = shift;
    my $self = { %Template };
    bless $self, $type;
    $self->{'display'} = [];
    $self->set_parameters(@_);
    $self;
}

#Called with the desired format to be parsed: either a string 'X' naming a 
#Parse::Format::X or a hint which will be recognised by that class.
sub parse {
    my ($self, $file, $format) = (shift, shift, shift);
    my ($library, $tmp, $bld, $aln, $dis, $header1, $header2, $header3, $loop);

    $self->set_parameters(@_);

    #load a parser for the desired format
    $tmp = "Bio::MView::Build::Format::$format";
    ($library = $tmp) =~ s/::/\//g;
    require "$library.pm";

    return  if $self->{'noparse'};    #load the parser, but do nothing

    $self->{'file'}   = $file;
    $self->{'format'} = lc $format;
    $self->{'class'}  = $tmp;

    #warn $self->{'format'}, "\n";

    no strict 'refs';
    $tmp = &{"${tmp}::parser"}();
    use strict 'refs';

    $self->{'stream'} = new NPB::Parse::Stream($file, $tmp);

    return undef  unless defined $self->{'stream'};

    ($loop, $header1, $header2, $header3) = (0, '', '', '');
    
    #$header1 = $self->header($self->{'quiet'});

    while (defined ($bld = $self->next)) {

        $bld->set_parameters(%{$self->{'bp'}});

        while (defined ($aln = $bld->next)) {

            next  unless $aln;    #null alignment
	    
	    if ($self->{'bp'}->{'mode'} eq 'rdb')
		{ print $bld->rdb($aln);     next }
	    if ($self->{'bp'}->{'mode'} eq 'pir')
		{ print $bld->pir($aln);     next }
	    if ($self->{'bp'}->{'mode'} eq 'msf')
		{ print $bld->msf($aln);     next }
	    if ($self->{'bp'}->{'mode'} eq 'pearson')
		{ print $bld->pearson($aln); next }
	    if ($self->{'bp'}->{'mode'} eq 'plain')
		{ print $bld->plain($aln);   next }
    
	    $dis = $self->add_display($bld, $aln);

	    if ($loop++ < 1) {
		$header2 = $bld->header($self->{'quiet'}) . $aln->header($self->{'quiet'});
	    }
	    $header3 = $bld->subheader($self->{'quiet'});

	    #add to display list
	    push @{$self->{'display'}}, [ $dis, $header1, $header2, $header3 ];

	    #display item now?
	    unless ($self->{'register'}) {
		$self->print;
		@{$self->{'display'}} = ();  #garbage collect
		#Universal::vmstat("print done (Manager)");
	    }

	    $header1 = $header2 = $header3 = '';

	    #drop old Align and Display objects: GC *before* next iteration!
	    $aln = $dis = undef;
        }

	#drop old Build object: GC *before* next iteration!
	$bld = undef;
    }
    $self;
}

sub set_parameters {
    my $self = shift;
    my ($key, $val);
    while (@_) {
        ($key, $val) = (shift, shift);
	#warn "($key, $val)\n";
        if (exists $Template{$key}) {
            $self->{$key} = $val;
            next;
        }
        if (exists $Template_Build_Param{$key}) {
            $self->{'bp'}->{$key} = $val;
            next;
        }
        if (exists $Template_Align_Param{$key}) {
            $self->{'ap'}->{$key} = $val;
            next;
        }
        warn "Bio::MView::Manager: unknown parameter '$key'\n";
    }
    #hack to share certain state
    if (exists $self->{'ap'}->{'label4'}) {
	$self->{'bp'}->{'showpcid'} = $self->{'ap'}->{'label4'};
    }
    $self;
}

#return next entry worth of parse data as in a Bio::MView::Build object 
#ready for parsing, or undef if no more data.
sub next {
    my $self = shift;
    my ($entry, $tmp);

    #free the last entry and garbage its Bio::MView::Build
    if (defined $self->{'filter'}) {
        $self->{'filter'}->get_entry->free;
        $self->{'filter'} = undef;
    }

    #read the next chunk of data
    $entry = $self->{'stream'}->get_entry;
    if (! defined $entry) {
        $self->{'stream'}->close;
        return undef;
    }

    #construct a new Bio::MView::Build
    return $self->{'filter'} = $self->{'class'}->new($entry);
}


#construct a header string describing this alignment
sub header {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s    if $quiet;
    $s .= "File: $self->{'file'}  Format: $self->{'format'}\n";
    Bio::MView::Display::displaytext($s);
}

sub add_display {
    my ($self, $bld, $aln) = @_;
    my ($dis, $tmp);

    my $ref = $bld->get_row_id($self->{'bp'}->{'ref_id'});

    #allow the Build instance to override the normal parameter
    #settings and to substitute specialised handlers for
    #'peculiar' alignments, eg., sequence versus secondary structure.
    $self->set_parameters($bld->change_parameters());
    $bld->change_alignment_type($aln);

    #Universal::vmstat("display constructor");
    $dis = new Bio::MView::Display($aln->init_display);
    #Universal::vmstat("display constructor DONE");
    
    #attach a ruler?
    if ($self->{'bp'}->{'ruler'}) {
        $tmp = $aln->build_ruler;
	$tmp->append_display($dis);
        #Universal::vmstat("ruler added");
    }

    #attach the alignment
    if ($self->{'bp'}->{'alignment'}) {
        $aln->set_color_scheme
	    (
	     'ref_id'      => $ref,
	     'coloring'    => $self->{'ap'}->{'aln_coloring'},
	     'colormap'    => $self->{'ap'}->{'aln_colormap'},
	     'colormap2'   => $self->{'ap'}->{'con_colormap'},
	     'group'       => $self->{'ap'}->{'aln_groupmap'},
	     'threshold'   => $self->{'ap'}->{'aln_threshold'},
	     'ignore'      => $self->{'ap'}->{'aln_ignore'},
	     'con_gaps'    => $self->{'ap'}->{'con_gaps'},
	     'css1'        => $self->{'ap'}->{'css1'},
	     'alncolor'    => $self->{'ap'}->{'alncolor'},
	     'labcolor'    => $self->{'ap'}->{'labcolor'},
	     'symcolor'    => $self->{'ap'}->{'symcolor'},
	     'gapcolor'    => $self->{'ap'}->{'gapcolor'},
	     'find'        => $self->{'ap'}->{'find'},
	    );
        #Universal::vmstat("set_color_scheme done");

	#determine GC policy by state of 'consensus': set means don't GC
	$aln->append_display($dis, $self->{'bp'}->{'consensus'});
        #Universal::vmstat("alignment added");
    }

    #attach consensus alignments?
    if ($self->{'bp'}->{'consensus'}) {
	$tmp = $aln->build_consensus_rows(
                                          $self->{'ap'}->{'con_groupmap'},
                                          $self->{'ap'}->{'con_threshold'},
                                          $self->{'ap'}->{'con_ignore'},
                                          $self->{'ap'}->{'con_gaps'},
                                         );
        $tmp->set_color_scheme(
			       'coloring'  => $self->{'ap'}->{'con_coloring'},
                               'colormap'  => $self->{'ap'}->{'aln_colormap'},
                               'colormap2' => $self->{'ap'}->{'con_colormap'},
                               'group'     => $self->{'ap'}->{'con_groupmap'},
                               'threshold' => $self->{'ap'}->{'con_threshold'},
                               'ignore'    => $self->{'ap'}->{'con_ignore'},
			       'css1'      => $self->{'ap'}->{'css1'},
                              );
	$tmp->append_display($dis);
        #Universal::vmstat("consensi added");
    }

    $dis;
}

#wrapper functions
sub check_identity_mode   { Bio::MView::Build::check_identity_mode(@_) }
sub check_display_mode 	  { Bio::MView::Build::check_display_mode(@_) }
sub check_hsp_tiling 	  { Bio::MView::Build::check_hsp_tiling(@_) }
sub check_molecule_type	  { Bio::MView::Align::check_molecule_type(@_) }
sub check_alignment_color_scheme { Bio::MView::Align::check_alignment_color_scheme(@_) }
sub check_consensus_color_scheme { Bio::MView::Align::check_consensus_color_scheme(@_) }
sub check_colormap     	  { Bio::MView::Align::check_colormap(@_) }
sub check_groupmap    	  { Bio::MView::Align::Consensus::check_groupmap(@_) }
sub check_ignore_class 	  { Bio::MView::Align::Consensus::check_ignore_class(@_) }
sub get_default_colormaps { Bio::MView::Align::get_default_colormaps(@_) }
sub get_default_groupmap  { Bio::MView::Align::Consensus::get_default_groupmap(@_) }
sub get_default_find_colormap
                          { Bio::MView::Align::get_default_find_colormap(@_) }
sub load_colormaps     	  { Bio::MView::Align::load_colormaps(@_) }
sub list_colormaps     	  { Bio::MView::Align::list_colormaps(@_) }
sub load_groupmaps     	  { Bio::MView::Align::Consensus::load_groupmaps(@_) }
sub list_groupmaps   	  { Bio::MView::Align::Consensus::list_groupmaps(@_) }
sub list_css              { Bio::MView::Align::print_css1_colormaps(@_) }

sub print {
    my ($self, $stm) = (@_, \*STDOUT);

    $self->{'posnwidth'} = 0;
    $self->{'labwidth0'} = 0;
    $self->{'labwidth1'} = 0;
    $self->{'labwidth2'} = 0;
    $self->{'labwidth3'} = 0;
    $self->{'labwidth4'} = 0;
    $self->{'labwidth5'} = 0;
    $self->{'labwidth6'} = 0;

    #consolidate field widths
    foreach (@{$self->{'display'}}) {
        $self->{'posnwidth'} = $_->[0]->{'posnwidth'}
            if $_->[0]->{'posnwidth'} > $self->{'posnwidth'};
        $self->{'labwidth0'} = $_->[0]->{'labwidth0'}
            if $_->[0]->{'labwidth0'} > $self->{'labwidth0'};
        $self->{'labwidth1'} = $_->[0]->{'labwidth1'}
            if $_->[0]->{'labwidth1'} > $self->{'labwidth1'};
        $self->{'labwidth2'} = $_->[0]->{'labwidth2'}
            if $_->[0]->{'labwidth2'} > $self->{'labwidth2'};
        $self->{'labwidth3'} = $_->[0]->{'labwidth3'}
            if $_->[0]->{'labwidth3'} > $self->{'labwidth3'};
        $self->{'labwidth4'} = $_->[0]->{'labwidth4'}
            if $_->[0]->{'labwidth4'} > $self->{'labwidth4'};
        $self->{'labwidth5'} = $_->[0]->{'labwidth5'}
            if $_->[0]->{'labwidth5'} > $self->{'labwidth5'};
        $self->{'labwidth6'} = $_->[0]->{'labwidth6'}
            if $_->[0]->{'labwidth6'} > $self->{'labwidth6'};
    }

    my $first = 1;
    #output
    while ($_ = shift @{$self->{'display'}}) {
	#Universal::vmstat("display");
	if ($self->{'html'}) {
	    my $s = '';
	    #body tag
	    if (! $self->{'ap'}->{'css1'}) {
                #supported in HTML 3.2:
		$s .= " BGCOLOR='$self->{'ap'}->{'alncolor'}'"
		    if defined $self->{'ap'}->{'alncolor'};
                #NOT supported in HTML 3.2? do them anyway:
		$s .= " TEXT='$self->{'ap'}->{'symcolor'}'"
		    if defined $self->{'ap'}->{'symcolor'};
		$s .= " LINK='$self->{'linkcolor'}'"
		    if defined $self->{'linkcolor'};
		$s .= " ALINK='$self->{'alinkcolor'}'"
		    if defined $self->{'alinkcolor'};
		$s .= " VLINK='$self->{'vlinkcolor'}'"
		    if defined $self->{'vlinkcolor'};
	    }
	    print $stm "<P>\n"  unless $first;
	    print $stm "<TABLE BORDER=0$s>\n";
	    #header
	    print $stm "<TR><TD><PRE>\n";
	    print $stm ($_->[1] ? $_->[1] : '');
	    print $stm ($_->[2] ? $_->[2] : '');
	    print $stm "</PRE></TD></TR>\n";
	    #subheader
	    if ($_->[3]) {
		print $stm "<TR><TD><PRE>\n";
		print $stm $_->[3];
		print $stm "</PRE></TD></TR>\n";
	    }
	    #alignment start
	    print $stm "<TR><TD>\n";
	} else {
	    #header
	    print $stm "\n"           if $_->[1] or $_->[2];
	    print $stm $_->[1],       if $_->[1];
	    print $stm $_->[2]        if $_->[2];
	    print "\n";
	    print $stm $_->[3], "\n"  if $_->[3];
	}
	#alignment
	$_->[0]->display(
			 'stream'    => $stm,
			 'html'      => $self->{'html'},
			 'bold'      => $self->{'ap'}->{'bold'},
			 'col'       => $self->{'ap'}->{'width'},
			 'label0'    => $self->{'ap'}->{'label0'},
			 'label1'    => $self->{'ap'}->{'label1'},
			 'label2'    => $self->{'ap'}->{'label2'},
			 'label3'    => $self->{'ap'}->{'label3'},
			 'label4'    => $self->{'ap'}->{'label4'},
			 'label5'    => $self->{'ap'}->{'label5'},
			 'label6'    => $self->{'ap'}->{'label6'},
			 'posnwidth' => $self->{'posnwidth'},
			 'labwidth0' => $self->{'labwidth0'},
			 'labwidth1' => $self->{'labwidth1'},
			 'labwidth2' => $self->{'labwidth2'},
			 'labwidth3' => $self->{'labwidth3'},
			 'labwidth4' => $self->{'labwidth4'},
			 'labwidth5' => $self->{'labwidth5'},
			 'labwidth6' => $self->{'labwidth6'},
			);
	if ($self->{'html'}) {
	    #alignment end
	    print $stm "</TD></TR>\n";
	    print $stm "</TABLE>\n";
	}
	#Universal::vmstat("display done");
	$_->[0]->free;
	#Universal::vmstat("display free done");

	$first = 0;
    }
    $self;
}


###########################################################################
1;
