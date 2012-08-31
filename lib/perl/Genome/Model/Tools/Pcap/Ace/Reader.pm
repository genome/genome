# 
# Copyright (C) 2004 Washington University Genome Sequencing Center
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#

# Last modified: <Wed, 2006/01/04 13:55:25 lcarmich linus58>

package Genome::Model::Tools::Pcap::Ace::Reader;
our $VERSION = 0.01;

=pod

=head1 NAME

AceReader - Ace file iterator

=head1 SYNOPSIS

    my $reader = new Genome::Model::Tools::Pcap::Ace::Reader(input => \*STDIN);
    while (my $obj = $reader->next_object()) {
        if ($obj->{'type'} eq 'contig') {
            ...
        }
        ...
    }

=head1 DESCRIPTION

Genome::Model::Tools::Pcap::Ace::Reader iterates over an ace file, returning one element at a time.

=head1 METHODS

=cut

use strict;
use warnings;
use Carp;

my $pkg = 'Genome::Model::Tools::Pcap::Ace::Reader';

=pod

=item new 

    my $reader = new Genome::Model::Tools::Pcap::Ace::Reader(\*STDIN);

=cut
sub new {
    croak("$pkg:new:no class given, quitting") if @_ < 1;
    my ($caller, $input) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = {};
    bless ($self, $class);
    
    %{$self->{'object_builders'}} = ( 
        AS => \&_build_assembly,
        AF => \&_build_read_position,
        CO => \&_build_contig,
        RD => \&_build_read,
        WA => \&_build_assembly_tag,
        CT => \&_build_contig_tag,
        RT => \&_build_read_tag,
        BS => \&_build_base_segment,
    );

    $self->{'input'} = $input;

    return $self;
}

=pod

=item next_object 

    my $obj_hashref = $reader->next_object();

    $obj_hashref->{'type'} eq 'contig'

    next_object returns the next object found in the ace file.  The return value is a
    hashref containing a 'type' key, and various other keys depending on the type:

    type eq 'assembly'
        contig_count
        read_count

    type eq 'contig'
        name
        base_count
        read_count
        base_seg_count
        u_or_c
        consensus
        base_qualities

    type eq 'read_position'
        read_name
        u_or_c
        position

    type eq 'base_segment'
        start_pos
        end_pos
        read_name

    type eq 'read'
        name
        padded_base_count
        info_count
        tag_count
        sequence
        qual_clip_start
        qual_clip_end
        align_clip_start
        align_clip_end
        description        - A hashref containing details about the trace
            CHROMAT_FILE
            PHD_FILE
            TIME

    type eq 'assembly_tag'
        tag_type
        program
        date
        data

    type eq 'contig_tag'
        contig_name
        tag_type
        program
        start_pos
        end_pos
        date
        no_trans
        data

    type eq 'read_tag'
        read_name
        tag_type
        program
        start_pos
        end_pos    
        date
        data
        

=cut

sub input_handle {
    my $self = shift;
    return $self->{input};
}

sub next_object {
    my ($self) = @_;
    my $IN = $self->{'input'};
    my $ret_val;
    while (my $line = <$IN>) {
        if ($line =~ /^\s*$/) {
            next;
        }
        chomp $line;
        my @tokens = split(/[ {]/,$line);
        if (@tokens > 0) {
            my $type = shift @tokens;
            if (exists $self->{'object_builders'}->{$type}) {
                $ret_val = $self->{'object_builders'}->{$type}->($self,$IN,\@tokens);
                $self->_fire_object_callback($ret_val);
                return $ret_val;
            }
        }
    }
    return undef;
}

sub parse {
    my $self = shift;
    while (my $obj = $self->next_object) {}
}

sub add_object_callback {
    my $self = shift;
    my $callback = shift;
    push @{$self->{callbacks}}, $callback;
}

sub _fire_object_callback {
    my $self = shift;
    my $object = shift;
    foreach my $callback (@{$self->{callbacks}}) {
        $callback->($object);
    }
}

sub width {
    my ($self,$width) = @_;
    if (defined $width) {
        $self->{'width'} = $width;
    }
    return $self->{'width'};
}

sub _build_assembly {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val = (
        type => 'assembly',
        contig_count => $token_ary_ref->[0],
        read_count => $token_ary_ref->[1],
    );
    return \%ret_val;
}

sub _build_contig {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val = (
        type => 'contig',
        name => $token_ary_ref->[0],
        base_count => $token_ary_ref->[1],
        read_count => $token_ary_ref->[2],
        base_seg_count => $token_ary_ref->[3],
        u_or_c => $token_ary_ref->[4],
    );

    my $consensus;

    my $line;
    while ($line = <$IN>) {
        if ($line =~ /^\s*$/) {
            last;
        }
        chomp $line;
        if (!$self->width) {
            $self->width(length $line);
        }
        $consensus .= $line;
    }
    $ret_val{'consensus'} = $consensus;
    while ($line = <$IN>) {
        if ($line =~ /^BQ/) {
            last;
        }
    }
    my @bq;
    while ($line = <$IN>) {
        if ($line =~ /^\s*$/) {
            last;
        }
        chomp $line;
        $line =~ s/^ //; # get rid of leading space
        push @bq, split(/ /,$line);
    }
    $ret_val{'base_qualities'} = \@bq;
    return \%ret_val;
}

sub _build_read_position {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val = (
        type => 'read_position',
        read_name => $token_ary_ref->[0],
        u_or_c => $token_ary_ref->[1],
        position => $token_ary_ref->[2],
    );
    return \%ret_val;
}

sub _build_base_segment {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val = (
        type => 'base_segment',
        start_pos => $token_ary_ref->[0],
        end_pos => $token_ary_ref->[1],
        read_name => $token_ary_ref->[2],
    );
    return \%ret_val;
}

sub _build_read {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val = (
        type => 'read',
        name => $token_ary_ref->[0],
        padded_base_count => $token_ary_ref->[1],
        info_count => $token_ary_ref->[2],
        tag_count => $token_ary_ref->[3],
    );
    my $sequence;
    my $line;
    while ($line = <$IN>) {
        if ($line =~ /^\s*$/) {
            last;
        }
        chomp $line;
        $sequence .= $line;
    }
    #my ($seq, $pads) = $self->un_pad_sequence($sequence);
    $ret_val{'sequence'} = $sequence;
    #$ret_val{'pads'} = $pads;
    while ($line = <$IN>) {
        chomp $line;
        if ($line =~ /^QA/) { my @tokens = split(/ /,$line); $ret_val{'qual_clip_start'} = $tokens[1]; $ret_val{'qual_clip_end'} = $tokens[2];
            $ret_val{'align_clip_start'} = $tokens[3];
            $ret_val{'align_clip_end'} = $tokens[4];
        }
        elsif ($line =~ /^DS/) {
            $line =~ s/ (\w+): /|$1|/g; #delimit the key-value pairs
            my @tokens = split(/\|/, $line); 
            shift @tokens; # drop the DS tag
            my %description = @tokens;
            $ret_val{'description'} = \%description;
            last;
        }
    }
    return \%ret_val;
}

sub _build_assembly_tag {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val;
    $ret_val{'type'} = 'assembly_tag';
    my $line = <$IN>;
    chomp $line;
	$line =~ s/^\s*// if $line =~ /\w/;;
    @ret_val{'tag_type', 'program', 'date'} = split(/ /, $line);
    my $data;
    while ($line = <$IN>) {
		$line =~ s/^\s*// if $line =~ /\w/;;
        if ($line =~ /^}/) {
            last;
        }
        $data .= $line;
    }
    $ret_val{'data'} = $data;

    return \%ret_val;
}

sub _build_contig_tag {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val;
    $ret_val{'type'} = 'contig_tag';
    my $line = <$IN>;
    chomp $line;
	$line =~ s/^\s*// if $line =~ /\w/;     
    @ret_val{'contig_name', 'tag_type', 'program', 'start_pos', 'end_pos', 'date', 'no_trans'} = split(/ /, $line);
    my $data;
    while ($line = <$IN>) {
		$line =~ s/^\s*// if $line =~ /\w/;
        if ($line =~ /^}/) {
            last;
        }
        $data .= $line;
    }
    $ret_val{'data'} = $data;

    return \%ret_val;
}

sub _build_read_tag {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val;
    $ret_val{'type'} = 'read_tag';
    my $line = <$IN>;
    chomp $line;
	$line =~ s/^\s*// if $line =~ /\w/;;
    @ret_val{'read_name', 'tag_type', 'program', 'start_pos', 'end_pos', 'date'} = split(/ /, $line);
    
    while (my $nextline= <$IN>){
	last if $nextline=~/^\s*}\s*\n?$/;
	$ret_val{data}.=$nextline;
    }
    return \%ret_val;
}

1;
#$Header$
