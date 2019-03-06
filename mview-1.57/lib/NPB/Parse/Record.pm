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

# $Id: Record.pm,v 1.55 2013/09/09 21:31:04 npb Exp $

###########################################################################
package NPB::Parse::Record;

use vars qw(@ISA $PackDelim $KeysDelim);
use NPB::Parse::Message;
use strict;

@ISA = qw(NPB::Parse::Message);

my $ignore_attributes = "text|offset|bytes|record_by_posn|record_by_type|indices|relative_key|absolute_key";

$PackDelim = "\0";
$KeysDelim  = "::";

# Warning! the {'record_by_*'} fields seem to require explicit dereferencing
# otherwise perl won't garbaged collect them until after the program exits, 
# even though they aren't circular. 
# Therefore, the top-level Record of a Record tree should call the free()
# method when the caller has finished with it a record, to allow normal
# garbage collection at run-time.
sub new {
    shift;                 #discard system supplied type
    my ($type, $parent, $text, $offset, $bytes) = (@_, -1, -1);

    my $self = {};
    bless $self, $type;    #use type specified by caller

    #warn $text->substr(), "\n";

    $self->{'text'}           = $text;
    $self->{'offset'}         = $offset;    #my string offset
    $self->{'bytes'}          = $bytes;     #my string length
    $self->{'record_by_posn'} = [];         #list of records
    $self->{'record_by_type'} = {};         #hash of record types
    $self->{'indices'}        = [];         #additional keys for indexing
    
    #subrecord counter
    $self->{'index'}         = $self->get_record_number($parent);

    #relative key for reporting
    $self->{'relative_key'}  = $self->get_type . $KeysDelim . $self->{'index'};

    #absolute hierarchical key for indexing/reporting
    $self->{'absolute_key'}  = '';
    $self->{'absolute_key'}  = $parent->{'absolute_key'}. $KeysDelim
	if defined $parent;
    $self->{'absolute_key'} .= $self->{'relative_key'};
    
    $self;
}

#sub DESTROY { warn "DESTROY $_[0]\n" } 

#Explicitly call this when finished with a Record to ensure indirect
#circular references through 'parent' and 'record_by_*' fields are
#broken by resetting the 'parent' field. Likewise, set the shared 'text'
#field to undefined. Without this, the Record hierarchy is never garbaged until
#after the program exits! The caller can supply a list of subrecord types, in
#which case only those will be marked for destruction, not the parent itslf.
sub free {
    my $self = shift;
    my ($type, $rec);

    #warn "FREE $self (@_)\n";

    if (@_) {
	#only free these subrecord types
	foreach $type (@_) {
	    if (exists $self->{'record_by_type'}->{$type}) {
		foreach $rec (@{$self->{'record_by_type'}->{$type}}) {
		    if (@$rec > 3) {
			$rec->[3]->free;    #recursively free subrecords
			splice @$rec, 3;    #remove parse object
		    }
		}
	    }
	}
	return $self;
    }

    #free every subrecord
    foreach $rec (@{$self->{'record_by_posn'}}) {
	if (@$rec > 3) {
	    #warn "  ", join(", ", @$rec), "\n";
	    $rec->[3]->free;                #recursively free subrecords
	    splice @$rec, 3;                #remove parse object
	}
    }
    
    #free subrecords
    $self->{'record_by_type'} = undef;
    $self->{'record_by_posn'} = undef;
    $self->{'text'} = undef;

    return $self;
}

#Given a triple (key, offset, bytecount), store these as an anonymous array
#under the appropriate keys in the record_by_ attributes.
sub add_record {
    my $self = shift;
    my $rec = [ @_ ];
    push @{$self->{'record_by_posn'}},              $rec;
    push @{$self->{'record_by_type'}->{$rec->[0]}}, $rec;
}

#Given a triple (key, offset, bytecount) in $rec, extract the record
#string, generate an object subclassed from $class, store this in $rec
#after the existing triple, and return the object.
sub add_object {
    my ($self, $rec) = @_;
    my ($class, $ob);
    $class = ref($self) . '::' . $rec->[0];
    #warn "ADD_OBJECT($class)\n";
    $ob = $class->new($self, $self->{'text'}, $rec->[1], $rec->[2]);
    push @$rec, $ob;
    $ob;
}

sub get_offset   { $_[0]->{'offset'} }
sub get_bytes    { $_[0]->{'bytes'} }
sub get_pos      { ($_[0]->{'offset'}, $_[0]->{'bytes'}) }
sub get_class    { ref($_[0]) }
sub get_index    { $_[0]->{'index'} }
sub relative_key { $_[0]->{'relative_key'} }
sub absolute_key { $_[0]->{'absolute_key'} }

sub get_type {
    my @id = split('::', ref($_[0]));
    pop @id;
}

sub get_record_number {
    my ($self, $parent) = (@_, undef);
    my ($i, $type);
    if (defined $parent) {
	$type   = $self->get_type;
	for ($i=0; $i < @{$parent->{'record_by_type'}->{$type}}; $i++) {
	    if ($parent->{'record_by_type'}->{$type}->[$i]->[1] == 
		$self->{'offset'}) {
		return $i+1;
	    }
	}
    }
    #there is no record number
    return 0;
}

sub push_indices {
    my ($self, @indices) = @_;
    push @{$self->{'indices'}}, @indices;
    $self;
}

#Return list of indices.
sub get_indices {
    @{$_[0]->{'indices'}}
}

#Return list of attributes of Record, excepting housekeeping ones.
sub list_attrs {
    my $self = shift;
    my ($key, @attr);
    foreach $key (grep !/^($ignore_attributes)$/o , keys %$self) {
	push @attr, $key;
    }
    @attr;
}

#Return list of record types in Record
sub list_record_types {
    sort keys %{$_[0]->{'record_by_type'}};
}

sub test_args {
    my $self = shift;
    my $line = shift;
    my $i; foreach $i (@_) {
	next    if $i ne '';
	$self->warn("incomplete match of line: $line");
    }
}

sub test_records {
    my $self = shift;
    local $_;
    foreach (@_) {
	unless (exists $self->{'record_by_type'}->{$_}) {
	    $self->warn("corrupt or missing '$_' record\n")
	}
    }
    $self;
}

sub pack_hash {
    my $self = shift;
    my ($val, $rec, @indicies);

    @indicies = $self->get_indices;

    $val = join($PackDelim,
		(
		 $self->relative_key,
		 $self->{'offset'},
		 $self->{'bytes'},
		 scalar(@indicies),
		 @indicies,
		));
    foreach $rec (@{$self->{'record_by_posn'}}) {
	#append (key, offset, bytes)
	$val .= join($PackDelim,
		     (
		      '',
		      $rec->[0],
		      $rec->[1],
		      $rec->[2],
		     ));
    }
    $val;
}

sub unpack_hash {
    my ($self, $val) = @_;
    my (@tmp, $num, $rec);

    @tmp = split($PackDelim, $val);

    $self->{'relative_key'} = shift @tmp;
    $self->{'offset'}       = shift @tmp;
    $self->{'bytes'}        = shift @tmp;
    $num                    = shift @tmp;
    push @{$self->{'indices'}}, splice(@tmp, 0, $num);
    while (@tmp) {
	$self->add_record(splice(@tmp, 0, 3));
    }
}

sub print {
    my ($self, $indent) = (@_, 0);
    my ($tmp, $r, $i, $rec) = ('');
    my $x = ' ' x $indent;

    printf "%sKey:   %s   Indices: [%s]\n", $x,
    $self->relative_key, join(',', $self->get_indices);

    printf "%sClass: %s\n", $x, $self;

    #print records in order of appearance in parent Record
    printf "%sSubrecords by posn:\n", $x;
    $self->print_records_by_posn($indent);

#    #print records in order of type
#    printf "%sSubrecords by type:\n", $x;
#    $self->print_records_by_type($indent);

    printf "%sMiscellaneous:\n", $x;
    printf "$x%20s -> %d\n",   'index', $self->{'index'};
    printf "$x%20s -> [%s]\n", 'pos', join(', ', $self->get_pos);
    if (defined $self->{'text'}) {
	## $tmp   = substr(${$self->{'text'}}, $self->{'offset'}, 30);
	$tmp   = $self->{'text'}->substr($self->{'offset'}, 30);
        ($tmp) = split("\n", $tmp);
    }
    printf "$x%20s -> \"%s\" ...\n",  'text', $tmp;
}

sub print_records_by_posn {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my ($rec, %count);

    foreach $rec (@{$self->{'record_by_posn'}}) {

	if (@{$self->{'record_by_type'}->{$rec->[0]}} > 1) {
	    printf "$x%20s -> [%s]\n",
	    $rec->[0] . '/' . ++$count{$rec->[0]}, join(", ", @$rec);
	} else {
	    printf "$x%20s -> [%s]\n",
	    $rec->[0],                             join(", ", @$rec);
	}
    }
}

sub print_records_by_type {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my ($key, $rec, $i);

    foreach $key (sort keys %{$self->{'record_by_type'}}) {
	if (@{$self->{'record_by_type'}->{$key}} > 1) {
	    for ($i=0; $i < @{$self->{'record_by_type'}->{$key}}; $i++) {
		$rec = $self->{'record_by_type'}->{$key}->[$i];
		printf "$x%20s -> [%s]\n",
		$rec->[0] . '/' . ($i+1),          join(", ", @$rec);
	    }
 	} else {
	    $rec = $self->{'record_by_type'}->{$key}->[0];
	    printf "$x%20s -> [%s]\n",
	    $rec->[0],                             join(", ", @$rec);
	} 
    }
}

#extract a substring with offset and bytecount from the 'text' attribute,
#defaulting to the entire string
sub substr {
    my $self = shift;
    my ($offset, $bytes) = (@_, $self->{'offset'}, $self->{'bytes'});

    #warn "Record::substr(@_)\n";

    if (! defined $self->{'text'}) {
	$self->die("substr() $self text field undefined");
    }

    ## substr(${$self->{'text'}}, $offset, $bytes);
    $self->{'text'}->substr($offset, $bytes);
}

#Given a record key for this entry, return an array of corresponding record
#strings or the first string if called in a scalar context. If called with
#just '*' return an array of all record strings at this level. If called 
#with no argument return the whole database entry string as a scalar.
sub string {
    my ($self, $key) = (@_, undef);
    my (@list, @keys, $rec) = ();

    #warn "string($key)\n";

    if (! defined $key) {
	return $self->substr($self->{'offset'}, $self->{'bytes'});
    }

    if ($key eq '*') {
	@keys = keys %{$self->{'record_by_type'}};
    }

    foreach $key (@keys) {
	foreach $rec ($self->key_range($key)) {
	    if (defined $rec->[3]) {
		#print "STRING() calling child\n";
		push @list, $rec->[3]->string;
	    } else {
		#print "STRING() calling substr [$rec->[1], $rec->[2]]\n";
		push @list, $self->substr($rec->[1], $rec->[2]);
	    }
	}
    }

    return @list    if wantarray;
    return $list[0] if @list;
    return '';  #no data
}

#Given a record key for this entry, parse corresponding records and return
#an array of the parsed objects or the first object if called in a scalar
#context. If called with no argument or just '*' parse everything at this 
#level.
sub parse {
    my ($self, $key) = (@_, '*');
    my (@list, @keys, $rec) = ();
    
    #warn "parse($key)\n";

    if ($key eq '*') {
	@keys = keys %{$self->{'record_by_type'}};
    } else {
	push @keys, $key;
    }
    
    foreach $key (@keys) {
	foreach $rec ($self->key_range($key)) {
	    if (defined $rec->[3]) {
		#object already processed
		push @list, $rec->[3];
	    } else {
		#build object and return for caller
		push @list, $self->add_object($rec);
	    }
	}
    }

    return @list    if wantarray;
    return $list[0] if @list;
    return undef;  #no data
}

#Split a key string of forms 'key' or 'key[N]' or 'key[M..N]' and return
#an array of records so indexed. Supplied range indices are assumed to be
#positive 1-based and are converted to 0-based for internal use.
sub key_range {
    my ($self, $key) = @_;
    my ($lo, $hi, @rec) = (undef, undef);

    if ($key =~ /(\S+)\s*\[\s*(\d+)\s*\]/) {
	($key, $lo, $hi) = ($1, $2, $2);
    }
    if ($key =~ /(\S+)\s*\[\s*(\d+)\s*\.\.\s*(\d+)\s*\]/) {
	if ($2 < $3) {
	    ($key, $lo, $hi) = ($1, $2, $3);
	} else {
	    ($key, $lo, $hi) = ($1, $3, $2);
	}
    }

    #does the key map to a record type?
    return unless exists $self->{'record_by_type'}->{$key};
	    
    #was a range requested?
    if (defined $lo) {

	$lo--; $hi--;    #convert to internal 0-based indices

	#want range?
	if ($lo != $hi) {
	    #adjust range to fit available records?
	    $lo = 0 
		unless defined $self->{'record_by_type'}->{$key}[$lo];
	    $hi = $#{$self->{'record_by_type'}->{$key}} 
	        unless defined $self->{'record_by_type'}->{$key}[$hi];
        } else {
	    #ignore non-existent single index
	    return unless $lo > -1 and defined $self->{'record_by_type'}->{$key}[$lo];
	}
	#get the slice
	@rec = @{$self->{'record_by_type'}->{$key}}[$lo..$hi];
    } else {
	#just get'em all
	@rec = @{$self->{'record_by_type'}->{$key}};
    }

    @rec;
}

#Given a record key for this entry, count how many corresponding records
#exist and return an array of counts or the first count if called in a
#scalar context. If called with no argument or just '*' count everything.
sub count {
    my ($self, $key) = (@_, '*');
    my (@list, @keys) = ();
    
    if ($key eq '*') {
	@keys = sort keys %{$self->{'record_by_type'}};
    } else {
	push @keys, $key;
    }
    
    foreach $key (@keys) {
	
	#is there a record instance for this type?
	if (exists $self->{'record_by_type'}->{$key}) {
	    push @list, scalar @{$self->{'record_by_type'}->{$key}};
	} else {
	    push @list, 0;
	}
    }
    
    return @list    if wantarray;
    return $list[0] if @list;
    return 0;  #no data
}

#tidy up identifiers: strip leading "/:" or ">" substrings
sub strip_leading_identifier_chars {
    my $string = shift;
    $string =~ s/^(\/:|>)//;
    $string;
}

#Returns $text less newlines while attempting to assemble hyphenated words
#and excess white space split over multiple lines correctly.
sub strip_english_newlines {
    my $text = shift;

    #multiple hyphens look like 'SO4--' or long dashes - word break with space
    $text =~ s/(--)\s*\n+\s*/$1 /sg;

    #single hyphens look like hyphenated words - join at the hyphen
    $text =~ s/(-)\s*\n+\s*/$1/sg;

    #remaining newlines - word break with space
    $text =~ s/\s*\n+\s*/ /sg;

    #skip trailing white space added after last newline was removed
    $text =~ s/\s+$//s;
    
    $text;
}

sub strip_leading_space {
    my $text = shift;
    $text =~ s/^[ \t]+//;
    $text;
}

sub strip_trailing_space {
    my $text = shift;
    $text =~ s/[ \t]+$//;
    $text;
}

sub strip_trailing_junk {
    my $text = shift;
    $text =~ s/[;:|,.-]+$//;  #default trailing PDB chain is often '_'
    $text;
}

sub clean_identifier {
    my $text = shift;
    $text = strip_leading_identifier_chars($text);
    $text = strip_trailing_junk($text);
    return $text;
}

#warn with error string: override Message::warn() to a) be quiet
#about the classname and b) interpolate NPB::Parse::Record->get_absolute_key.
sub warn {
    my $self = shift;
    my $s = $_[$#_]; chomp $s;
    my $t = $self; $t =~ s/=.*//;
    warn "Warning $t ($self->{'absolute_key'}) $s\n";
}


###########################################################################
package NPB::Parse::Record_Stream;

use vars qw(@ISA);
use NPB::Parse::Message;
use NPB::Parse::Substring;
use strict;

@ISA = qw(NPB::Parse::Message);

sub new {
    my $type = shift;
    my $self = {};
    my ($entry, $depth,  $text, $offset, $bytes) = (@_, 0);
    $self->{'entry'}  = $entry;
    $self->{'depth'}  = $depth;
    if (defined $text and ref $text) {
	#use supplied text and positions (or defaults thereof)
	$self->{'text'}   = new NPB::Parse::Substring($text);
	$self->{'offset'} = defined $offset ? $offset : 0;
	$self->{'bytes'}  = defined $bytes ? $bytes : length($$text);
    } else {
	#use entry object's text and positions
	$self->{'text'}   = $entry->{'text'};
	$self->{'offset'} = $entry->{'offset'};
	$self->{'bytes'}  = $entry->{'bytes'};
    }
    bless $self, $type;
    $self->reset;
}

sub DESTROY { 
    #warn "DESTROY $_[0]\n";
    ####NPB no effect
    map { $_[0]->{$_} = undef } keys %{$_[0]};
}

sub reset {
    $_[0]->{'cursor'} = $_[0]->{'offset'};
    $_[0]->{'limit'}  = $_[0]->{'offset'} + $_[0]->{'bytes'};
    $_[0]->{'line'}   = '';
    $_[0]->{'length'} = 0;
    #warn "INITIAL=(c=$_[0]->{'cursor'} l=$_[0]->{'limit'})\n"; #NPB
    $_[0];
}

sub backup {
    my $self = shift;
    my $bytes = $self->{'depth'} + length($self->{'line'});
    $self->{'cursor'} -= $self->{'length'};
    $self->{'line'}    = '';
    $self->{'length'}  = 0;
}

#return the remaining unprocessed stream without altering the cursor
sub inspect_stream {
    return ''    if $_[0]->{'cursor'} >= $_[0]->{'limit'};
    $_[0]->{'text'}->substr($_[0]->{'cursor'},
			    $_[0]->{'limit'} - $_[0]->{'cursor'} + 1);
}

#return next line of text or undef if the text stream is done. chomp returned
#string if optional non-zero argument (defaults to zero).
sub next_line {
    my ($self, $chomp) = (@_, 0);

    if ($self->{'cursor'} >= $self->{'limit'}) {
	return $self->{'line'} = undef;    #so next call won't see it
    }

    #ignore 'depth' leading characters
    my $ptr = $self->{'cursor'} + $self->{'depth'};
    
    #how many bytes to read? upto newline or until limiting bytecount
    $self->{'line'} = $self->{'text'}->getline($ptr);
    my $bytes = length $self->{'line'};
    if ($self->{'cursor'} + $bytes > $self->{'limit'}) {
	$bytes = $self->{'limit'} - $ptr;
	$self->{'line'} = substr($self->{'line'}, 0, $bytes);
    }

    #how many bytes were actually read?
    $self->{'length'}  = $self->{'depth'} + $bytes;
    $self->{'cursor'} += $self->{'length'};

    #remove DOS CR; this should not be done at this level
    $self->{'line'} =~ s/\015//g;

    #return last thing read
    chomp $self->{'line'}    if $chomp;
    $self->{'line'};
}

#Read $count lines or all lines until (EOF or end-of-string) if $count==0.
#Store final $record in 'entry' using $key if set (defaults to unset).
#Assumes next_line() called just previously.
sub scan_lines {
    my ($self, $count, $key) = (@_, 0);
    my ($offset, $bytes, $line, $record) = (0, 0, '', '');
    
    $record = $self->{'line'};
    $offset = $self->{'cursor'} - $self->{'length'};

    #warn "scan_lines() looking at ($count) lines\n";

    $count--    if $record;    #already seen one line

    if ($count > 0) {
	#scan $count lines
	while ($count-- > 0 and defined ($line = $self->next_line)) {
	    $record .= $line;
	}
    } else {
	#scan everything
	while (defined ($line = $self->next_line)) {
	    $record .= $line;
	}
    }
    #no $self->backup as we've read exactly the right amount
    $bytes = $self->{'cursor'} - $offset;

    $self->{'entry'}->add_record($key, $offset, $bytes)    if $key;
    $record;
}

#Read >= 1 record lines terminating on failure to match $pattern. Store
#final $record in 'entry' using $key if set (defaults to unset).  Assumes
#next_line() called just previously.
sub scan_while {
    my ($self, $pattern, $key) = (@_, 0);
    my ($offset, $bytes, $line, $record) = (0, 0, '', '');

    #warn "scan_while() looking at /$pattern/\n";

    $record = $self->{'line'};
    $offset = $self->{'cursor'} - $self->{'length'};

    while (defined ($line = $self->next_line)) {
	if ($line =~ /$pattern/) {
	    $record .= $line;
	    next;
	}
	$self->backup;
	last;
    }
    $bytes = $self->{'cursor'} - $offset;

    $self->{'entry'}->add_record($key, $offset, $bytes)    if $key;
    $record;
}

#Read >= 1 record lines terminating on matching $pattern. Store final
#$record in 'entry' using $key if set (defaults to unset).  Assumes
#next_line() called just previously. Consumed record EXCLUDES matched line.
sub scan_until {
    my ($self, $pattern, $key) = (@_, 0);
    my ($offset, $bytes, $line, $record) = (0, 0, '', '');

    #warn "scan_until() looking until /$pattern/\n";

    $record = $self->{'line'};
    $offset = $self->{'cursor'} - $self->{'length'};

    while (defined ($line = $self->next_line)) {
	if ($line =~ /$pattern/) {
	    $self->backup;
	    last;
	}
	$record .= $line;
	next;
    }
    $bytes = $self->{'cursor'} - $offset;

    $self->{'entry'}->add_record($key, $offset, $bytes)    if $key;
    $record;
}

#Read >= 1 record lines terminating on matching $pattern. Store final
#$record in 'entry' using $key if set (defaults to unset).  Assumes
#next_line() called just previously. Consumed record INCLUDES matched line.
sub scan_until_inclusive {
    my ($self, $pattern, $key) = (@_, 0);
    my ($offset, $bytes, $line, $record) = (0, 0, '', '');

    #warn "scan_until_inclusive() looking until /$pattern/\n";

    $record = $self->{'line'};
    $offset = $self->{'cursor'} - $self->{'length'};

    while (defined ($line = $self->next_line)) {
	$record .= $line;
	last    if $line =~ /$pattern/;
	next;
    }
    $bytes = $self->{'cursor'} - $offset;

    $self->{'entry'}->add_record($key, $offset, $bytes)    if $key;
    $record;
}

#Read >= 1 record lines terminating on matching $pattern. Store final
#$record in 'entry' using $key if set (defaults to unset).  Assumes
#next_line() called just previously. Consumed record EXCLUDES matched line.
#Skips initial $skipcount instances of $pattern.
sub scan_skipping_until {
    my ($self, $pattern, $skip, $key) = (@_, 0);
    my ($offset, $bytes, $line, $record) = (0, 0, '', '');

    #warn "scan_skipping_until() looking until /$pattern/\n";

    $record = $self->{'line'};
    $offset = $self->{'cursor'} - $self->{'length'};

    while (defined ($line = $self->next_line)) {
	if ($line =~ /$pattern/) {
	    if ($skip-- < 1) {
	        $self->backup;
	        last;
            }
	}
	$record .= $line;
	next;
    }
    $bytes = $self->{'cursor'} - $offset;

    $self->{'entry'}->add_record($key, $offset, $bytes)    if $key;
    $record;
}

#Read >= 1 record lines terminating on failure to match empty lines or
#initial blank space up to $nest characters. Store final $record in 'entry'
#using $key if set (defaults to unset).  Assumes next_line() called just
#previously.
sub scan_nest {
    my ($self, $nest, $key) = (@_, 0);
    my ($offset, $bytes, $line, $record) = (0, 0, '', '');

    #warn "scan_nest() looking at nest depth $nest\n";

    $record = $self->{'line'};
    $offset = $self->{'cursor'} - $self->{'length'};

    while (defined ($line = $self->next_line)) {
	if ($line =~ /^(\s{$nest}|$)/) {
	    $record .= $line;
	    next;
	}
	$self->backup;
	last;
    }
    $bytes = $self->{'cursor'} - $offset;

    $self->{'entry'}->add_record($key, $offset, $bytes)    if $key;
    $record;
}


###########################################################################
1;
