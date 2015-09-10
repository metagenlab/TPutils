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

# $Id: FASTA3X.pm,v 1.3 2013/10/20 22:01:50 npb Exp $

###########################################################################
#
# FASTA 3X (34/35/36)
#
#   fasta, fastx, fasty, tfasta, tfastx, tfasty, tfastxy,
#   fastm, fastf, fasts, ggsearch, glsearch, ssearch
#
###########################################################################
###########################################################################
use Bio::MView::Build::Format::FASTA3;


###########################################################################
package Bio::MView::Build::Format::FASTA3X;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3);


###########################################################################
###########################################################################
package Bio::MView::Build::Row::FASTA3X;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::fasta;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3::fasta);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::fastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3::fastx);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::fasty;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3::fasty);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::tfasta;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3::tfasta);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::tfastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3::tfastx);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::tfasty;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3::tfasty);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::tfastxy;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3::tfastxy);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::ssearch;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA);

sub new {
    my $type = shift;
    my ($num, $id, $desc, $score, $bits, $e, $query_orient, $sbjct_orient)
	= @_;
    my $self = new Bio::MView::Build::Row($num, $id, $desc);
    $self->{'score'} 	    = $score;
    $self->{'bits'}    	    = $bits;
    $self->{'e'}       	    = $e;
    $self->{'query_orient'} = $query_orient;
    $self->{'sbjct_orient'} = $sbjct_orient;
    bless $self, $type;
}

sub data  {
    return sprintf("%5s %7s %9s %2s %2s",
		   'S-W', 'bits', 'E-value', 'qy', 'ht')
	unless $_[0]->num;
    return sprintf("%5s %7s %9s %2s %2s",
		   $_[0]->{'score'}, $_[0]->{'bits'}, $_[0]->{'e'},
		   $_[0]->{'query_orient'}, $_[0]->{'sbjct_orient'});
}

sub rdb {
    my ($self, $mode) = (@_, 'data');
    my $s = Bio::MView::Build::Row::rdb($self, $mode);
    return join("\t", $s, $self->{'score'}, $self->{'bits'}, $self->{'e'},
		$self->{'query_orient'}, $self->{'sbjct_orient'})
	if $mode eq 'data';
    return join("\t", $s, 'S-W', 'bits', 'E-value',
		'query_orient', 'sbjct_orient')
	if $mode eq 'attr';
    return join("\t", $s, '5N', '7N', '9N', '2S', '2S')
	if $mode eq 'form';
    '';
}

sub range {
    my $self = shift;
    $self->SUPER::range($self->{'query_orient'});
}

sub assemble { my $self = shift; $self->assemble_fasta(@_) }


###########################################################################
package Bio::MView::Build::Row::FASTA3X::fastm;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA);

sub new {
    my $type = shift;
    my ($num, $id, $desc, $initn, $init1, $bits, $e, $sn, $sl,
	$query_orient, $sbjct_orient) = @_;
    my $self = new Bio::MView::Build::Row($num, $id, $desc);
    $self->{'initn'} 	    = $initn;
    $self->{'init1'} 	    = $init1;
    $self->{'bits'}    	    = $bits;
    $self->{'e'}       	    = $e;
    $self->{'sn'}      	    = $sn;
    $self->{'sl'}      	    = $sl;
    $self->{'query_orient'} = $query_orient;
    $self->{'sbjct_orient'} = $sbjct_orient;
    bless $self, $type;
}

sub data  {
    return sprintf("%5s %5s %7s %9s %3s %3s %2s %2s",
		   'initn', 'init1', 'bits', 'E-value', 'sn', 'sl', 'qy', 'ht')
	unless $_[0]->num;
    return sprintf("%5s %5s %7s %9s %3s %3s %2s %2s",
		   $_[0]->{'initn'}, $_[0]->{'init1'}, $_[0]->{'bits'},
		   $_[0]->{'e'}, $_[0]->{'sn'}, $_[0]->{'sl'},
		   $_[0]->{'query_orient'}, $_[0]->{'sbjct_orient'});
}

sub rdb {
    my ($self, $mode) = (@_, 'data');
    my $s = Bio::MView::Build::Row::rdb($self, $mode);
    return join("\t", $s, $self->{'initn'}, $self->{'init1'}, $self->{'bits'},
		$self->{'e'}, $self->{'sn'}, $self->{'sl'},
		$self->{'query_orient'}, $self->{'sbjct_orient'})
	if $mode eq 'data';
    return join("\t", $s, 'initn', 'init1', 'bits', 'E-value', 'sn', 'sl',
		'query_orient', 'sbjct_orient')
	if $mode eq 'attr';
    return join("\t", $s, '5N', '5N', '7N', '9N', '3N', '3N', '2S', '2S')
	if $mode eq 'form';
    '';
}

sub range {
    my $self = shift;
    $self->SUPER::range($self->{'query_orient'});
}

sub assemble { my $self = shift; $self->assemble_fasta(@_) }


###########################################################################
package Bio::MView::Build::Row::FASTA3X::ggsearch;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3X::ssearch);

sub data  {
    return sprintf("%5s %7s %9s %2s %2s",
		   'N-W', 'bits', 'E-value', 'qy', 'ht')
	unless $_[0]->num;
    return sprintf("%5s %7s %9s %2s %2s",
		   $_[0]->{'score'}, $_[0]->{'bits'}, $_[0]->{'e'},
		   $_[0]->{'query_orient'}, $_[0]->{'sbjct_orient'});
}

sub rdb {
    my ($self, $mode) = (@_, 'data');
    my $s = Bio::MView::Build::Row::rdb($self, $mode);
    return join("\t", $s, $self->{'score'}, $self->{'bits'}, $self->{'e'},
		$self->{'query_orient'}, $self->{'sbjct_orient'})
	if $mode eq 'data';
    return join("\t", $s, 'N-W', 'bits', 'E-value',
		'query_orient', 'sbjct_orient')
	if $mode eq 'attr';
    return join("\t", $s, '5N', '7N', '9N', '2S', '2S')
	if $mode eq 'form';
    '';
}


###########################################################################
package Bio::MView::Build::Row::FASTA3X::glsearch;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3X::ggsearch);


###########################################################################
###########################################################################
package Bio::MView::Build::Format::FASTA3X::fasta;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3::fasta);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::fastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3::fastx);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::fasty;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3::fasty);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::tfasta;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3::tfasta);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::tfastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3::tfastx);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::tfasty;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3::tfasty);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::tfastxy;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3::tfastxy);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::ssearch;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA); #note

sub subheader {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s    if $quiet;
    $s  = $self->SUPER::subheader($quiet);
    $s .= "Query orientation: " . $self->strand . "\n";
    $s;
}

sub parse {
    my $self = shift;
    return $self->parse_body('ssearch', @_);
}

sub parse_body {
    my ($self, $hint) = (shift, shift);
    my ($match, $sum, $aln, $query, $key);
    my ($rank, $use, %hit, @hit) = (0);

    #the actual Row subclass to build
    my $class = "Bio::MView::Build::Row::FASTA3X::$hint";

    #identify the query itself
    $match = $self->{'entry'}->parse(qw(HEADER));

    #all strands done?
    return  unless defined $self->schedule_by_strand;

    #identify the query itself
    $match = $self->{'entry'}->parse(qw(HEADER));

    if ($match->{'query'} ne '') {
	$query = $match->{'query'};
    } elsif ($match->{'queryfile'} =~ m,.*/([^\.]+)\.,) {
	$query = $1;
    } else {
	$query = 'Query';
    }

    #fasta run with no hits
    my $rankparse = $self->{'entry'}->parse(qw(RANK));
    return []  unless defined $rankparse;

    push @hit, new $class(
	'',
	$query,
	'',
	'',
	'',
	'',
	$self->strand,
	'',
	);
    
    #extract cumulative scores and identifiers from the ranking
    foreach $match (@{ $rankparse->{'hit'} }) {

	$rank++;

	#check row wanted, by num OR identifier OR row count limit OR score:
	#in ssearch rankings, 'score' seems the same as 'opt in the summaries
	#so use the same fasta use_row filter
	last  if ($use = $self->use_row($rank, $rank, $match->{'id'},
					$match->{'score'})
		 ) < 0;
	next  unless $use;

	#warn "KEEP: ($rank,$match->{'id'})\n";

	$key = $match->{'id'} . $match->{'score'} . $match->{'expect'};

	#warn "ADD: [$key]\n";

	push @hit, new $class(
	    $rank,
	    $match->{'id'},
	    $match->{'desc'},
	    $match->{'score'},
	    $match->{'bits'},
	    $match->{'expect'},
	    $self->strand,
	    '',
	    );
	$hit{$key} = $#hit;
    }

    #pull out each hit
    foreach $match ($self->{'entry'}->parse(qw(MATCH))) {

	#first the summary
	$sum = $match->parse(qw(SUM));

	#only read hits already seen in ranking
	while (1) {
	    #SSEARCH3X reports two s-w scores, either might match:
	    $key = $sum->{'id'} . $sum->{'opt'} . $sum->{'expect'};
	    last  if exists $hit{$key};
	    $key = $sum->{'id'} . $sum->{'score'} . $sum->{'expect'};
	    last  if exists $hit{$key};
	    $key = '';
	    last;
	}
	next  unless exists $hit{$key};
	#warn "SEE: [$key]\n";

	#override the row description
	if ($sum->{'desc'}) {
	    $hit[$hit{$key}]->{'desc'} = $sum->{'desc'};
	}

	#then the individual matched fragments
	foreach $aln ($match->parse(qw(ALN))) {

	    #ignore other query strand orientation
            next  unless $aln->{'query_orient'} eq $self->strand;

	    $aln = $match->parse(qw(ALN));
	    
	    #$aln->print;
	    
	    #for FASTA gapped alignments
	    $self->strip_query_gaps(\$aln->{'query'}, \$aln->{'sbjct'},
				    $aln->{'query_leader'},
                                    $aln->{'query_trailer'});
	    
	    $hit[0]->add_frag
		(
		 $aln->{'query'},
		 $aln->{'query_start'},
		 $aln->{'query_stop'},
		 $aln->{'query_start'},
		 $aln->{'query_stop'},
		 0,
		 0,
		);
	    
	    $hit[$hit{$key}]->add_frag
		(
		 $aln->{'sbjct'},
		 $aln->{'query_start'},
		 $aln->{'query_stop'},
		 $aln->{'query_start'},
		 $aln->{'query_stop'},
		 $aln->{'sbjct_start'},
		 $aln->{'sbjct_stop'},
		);

	    #override row data
	    $hit[$hit{$key}]->{'sbjct_orient'} = $aln->{'sbjct_orient'};
	}
    }

    $self->discard_empty_ranges(\@hit);

    #free objects
    $self->{'entry'}->free(qw(HEADER RANK MATCH));

    #map { $_->print; print "\n" } @hit;

    return \@hit;
}


###########################################################################
package Bio::MView::Build::Format::FASTA3X::ggsearch;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3X::ssearch); #note

sub parse {
    my $self = shift;
    return $self->parse_body('ggsearch', @_);
}


###########################################################################
package Bio::MView::Build::Format::FASTA3X::glsearch;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3X::ssearch); #note

sub parse {
    my $self = shift;
    return $self->parse_body('glsearch', @_);
}


###########################################################################
package Bio::MView::Build::Format::FASTA3X::fastm;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA); #note

sub subheader {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s    if $quiet;
    $s  = $self->SUPER::subheader($quiet);
    $s .= "Query orientation: " . $self->strand . "\n";
    $s;
}

sub parse {
    my $self = shift;
    return $self->parse_body('fastm', @_);
}

sub parse_body {
    my ($self, $hint) = (shift, shift);
    my ($match, $sum, $aln, $query, $key);
    my ($rank, $use, %hit, @hit) = (0);

    #the actual Row subclass to build
    my $class = "Bio::MView::Build::Row::FASTA3X::$hint";

    #identify the query itself
    $match = $self->{'entry'}->parse(qw(HEADER));

    #all strands done?
    return  unless defined $self->schedule_by_strand;

    #identify the query itself
    $match = $self->{'entry'}->parse(qw(HEADER));

    if ($match->{'query'} ne '') {
	$query = $match->{'query'};
    } elsif ($match->{'queryfile'} =~ m,.*/([^\.]+)\.,) {
	$query = $1;
    } else {
	$query = 'Query';
    }

    #fasta run with no hits
    my $rankparse = $self->{'entry'}->parse(qw(RANK));
    return []  unless defined $rankparse;

    push @hit, new $class(
	'',
	$query,
	'',
	'',
	'',
	'',
	'',
	'',
	'',
	$self->strand,
	'',
	);
    
    #extract cumulative scores and identifiers from the ranking
    foreach $match (@{ $rankparse->{'hit'} }) {

	$rank++;

	#check row wanted, by num OR identifier OR row count limit OR initn OR
	#initn in fastm rankings.
	last  if ($use = $self->use_row($rank, $rank, $match->{'id'},
					$match->{'initn'})
		 ) < 0;
	next  unless $use;

	#warn "KEEP: ($rank,$match->{'id'})\n";

	$key = $match->{'id'} . $match->{'initn'} . $match->{'expect'};

	#warn "ADD: [$key]\n";

	push @hit, new $class(
	    $rank,
	    $match->{'id'},
	    $match->{'desc'},
	    $match->{'initn'},
	    $match->{'init1'},
	    $match->{'bits'},
	    $match->{'expect'},
	    $match->{'sn'},
	    $match->{'sl'},
	    $self->strand,
	    '',
	    );
	$hit{$key} = $#hit;
    }

    #pull out each hit
    foreach $match ($self->{'entry'}->parse(qw(MATCH))) {

	#first the summary
	$sum = $match->parse(qw(SUM));

	#only read hits already seen in ranking
	while (1) {
	    #FASTM3X reports three s-w scores, any might match:
	    $key = $sum->{'id'} . $sum->{'opt'} . $sum->{'expect'};
	    last  if exists $hit{$key};
	    $key = $sum->{'id'} . $sum->{'initn'} . $sum->{'expect'};
	    last  if exists $hit{$key};
	    $key = $sum->{'id'} . $sum->{'init1'} . $sum->{'expect'};
	    last  if exists $hit{$key};
	    $key = '';
	    last;
	}
	next  unless exists $hit{$key};
	#warn "SEE: [$key]\n";

	#override the row description
	if ($sum->{'desc'}) {
	    $hit[$hit{$key}]->{'desc'} = $sum->{'desc'};
	}

	#then the individual matched fragments
	foreach $aln ($match->parse(qw(ALN))) {

	    #ignore other query strand orientation
            next  unless $aln->{'query_orient'} eq $self->strand;

	    $aln = $match->parse(qw(ALN));
	    
	    #$aln->print;
	    
	    #for FASTA gapped alignments
	    $self->strip_query_gaps(\$aln->{'query'}, \$aln->{'sbjct'},
				    $aln->{'query_leader'},
                                    $aln->{'query_trailer'});
	    
	    $hit[0]->add_frag
		(
		 $aln->{'query'},
		 $aln->{'query_start'},
		 $aln->{'query_stop'},
		 $aln->{'query_start'},
		 $aln->{'query_stop'},
		 0,
		 0,
		);
	    
	    $hit[$hit{$key}]->add_frag
		(
		 $aln->{'sbjct'},
		 $aln->{'query_start'},
		 $aln->{'query_stop'},
		 $aln->{'query_start'},
		 $aln->{'query_stop'},
		 $aln->{'sbjct_start'},
		 $aln->{'sbjct_stop'},
		);

	    #override row data
	    $hit[$hit{$key}]->{'sbjct_orient'} = $aln->{'sbjct_orient'};
	}
    }

    $self->discard_empty_ranges(\@hit);

    #free objects
    $self->{'entry'}->free(qw(HEADER RANK MATCH));

    #map { $_->print; print "\n" } @hit;

    return \@hit;
}


###########################################################################
package Bio::MView::Build::Format::FASTA3X::fasts;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3X::fastm); #note


###########################################################################
package Bio::MView::Build::Format::FASTA3X::fastf;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3X::fastm); #note


###########################################################################
1;
