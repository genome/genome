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

# Last modified: <Thu, 2008/01/31 16:11:56 ebelter linus108>

package Finishing::Assembly::Phd::Reader;
our $VERSION = 0.01;

=pod

=head1 NAME

PhdReader - Phd file iterator

=head1 SYNOPSIS

    my $reader = new Finishing::Assembly::Phd::Reader(\*STDIN);
    while (my $obj = $reader->nextObject()) {
    }

=head1 DESCRIPTION

Finishing::Assembly::Phd::Reader iterates over a phd file, returning one element at a time.

=head1 METHODS

=cut

use strict;
use warnings;
use Carp;

use Finishing::Assembly::Factory;

my $pkg = 'Finishing::Assembly::Phd::Reader';

=pod

=item new 

    my $reader = new Finishing::Assembly::Phd::Reader;

=cut
sub new {
    croak("$pkg:new:no class given, quitting") if @_ < 1;
    my ($caller, $arg) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = {};
    bless ($self, $class);
    
    %{$self->{'object_builders'}} = ( 
        BEGIN_SEQUENCE => \&read_sequence_name,
        BEGIN_COMMENT => \&read_comment,
        BEGIN_DNA => \&read_DNA,
        "WR{" => \&read_WR,
        BEGIN_TAG => \&read_tag,
    );

#    $self->{'input'} = $arg;

    return $self;
}

=pod

=item Finishing::Assembly::Phd::Reader::read 

    $phd = $reader->read(\*STDIN);

=cut
sub read {
    my ($self,$IN) = @_;
	my $factory = Finishing::Assembly::Factory->connect("source");
    my $phd = $factory->create_assembled_read();
    while (my $line = <$IN>) {
        chomp $line;
        my @tokens = split(/ /, $line);
        my $tag = shift @tokens;
        if (defined $tag && $tag && exists $self->{'object_builders'}->{$tag}) {
            $self->{'object_builders'}->{$tag}->($self,\@tokens, $IN, $phd);
        }
    }
    return $phd;
}

sub read_sequence_name {
    my ($self,$tokens, $IN, $phd) = @_;
    $phd->name($tokens->[0]);
}

sub read_comment {
    my ($self,$tokens, $IN, $phd) = @_;
    my %comments;
    while (my $line = <$IN>) {
        if ($line =~ /END_COMMENT/) {
            last;
        }
        if ($line !~ /^\W*$/) {
            chomp $line;
            my @element = split(/: /,$line);
            $comments{$element[0]} = $element[1];
        }
    }
    $phd->comments(\%comments);
}

sub read_DNA {
    my ($self,$tokens, $IN, $phd) = @_;
    my $bases;
    my @quality;
    my @positions;
    while (my $line = <$IN>) {
        if ($line =~ /END_DNA/) {
            last;
        }
        chomp $line;
        my @pos = split(/ /,$line);
        $bases .= $pos[0];
        push @quality, $pos[1];
        push @positions, $pos[2];
    }
    $phd->unpadded_base_string($bases);
    $phd->unpadded_base_quality(\@quality);
    $phd->unpadded_chromat_positions(\@positions);
}

sub read_WR {
    my ($self,$tokens, $IN, $phd) = @_;
    my @wr;
    if (defined $phd->wr) {
        @wr = @{$phd->wr};
    }
    my $wrText;
    while (my $line = <$IN> ) {
        last if $line =~ /^}/;
        $wrText .= $line;
    }
    push @wr, $wrText;
    $phd->wr(\@wr);
    $self->_parse_wr($phd, $wrText);
}

sub _parse_wr {
    my ($self, $phd, $wr) = @_;
    if ($wr =~ /template/s) {
        $wr =~ /name: (\S*)/s;
        $phd->template($1);
    }
    elsif ($wr =~ /primer/s) {
        $wr =~ /type: (\S*)/s;
        $phd->primer($1);
    }
}

sub read_tag {
    my ($self,$tokens, $IN, $phd) = @_;
    my %tag_data;
    my $text = undef;
    while (my $line = <$IN>) {
        chomp $line;
        last if $line =~ /^END_TAG/;
        if ($line =~ /BEGIN_COMMENT/) {
            while (my $comment_line = <$IN>) {
                last if $comment_line =~ /END_COMMENT/;
                $text .= $comment_line;
            }
        }
        my @element = split /: /, $line;
        $tag_data{$element[0]} = $element[1] if defined $element[0];
    }
    my @pos = split / /, $tag_data{UNPADDED_READ_POS};
    my $fact = Finishing::Assembly::Factory->connect("source");
	my $tag = $fact->create_read_tag(
        type   => $tag_data{TYPE},
        source => $tag_data{SOURCE},
        date   => $tag_data{DATE},
        start  => $pos[0],
        stop   => $pos[1],
        parent => $phd->name,
        scope  => 'Phd',
    );
    $tag->text($text) if defined $text;
    my @tags = @{$phd->tags};
    push @tags, $tag;
    $phd->tags(\@tags);
    #$phd->add_tag($tag);
}

1;
