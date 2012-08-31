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

# Last modified: <Wed, 2005/09/28 10:46:49 lcarmich linus58>

package Genome::Model::Tools::Pcap::Ace::Writer;
our $VERSION = 0.01;

=pod

=head1 NAME
Genome::Model::Tools::Pcap::Ace::Writer - Write an ace file one element at a time

=head1 SYNOPSIS

my $writer = new Genome::Model::Tools::Pcap::Ace::Writer(\*STDOUT);

foreach $contig_hashref (@contigs) {
    $writer->write_object($contig_hashref);
}

=head1 DESCRIPTION

Genome::Model::Tools::Pcap::Ace::Writer takes hashrefs representing elements of an assembly and writes them out
to an ace file.

=cut

use strict;
use warnings;
use Carp;
use Data::Dumper;

my $pkg = 'Genome::Model::Tools::Pcap::Ace::Writer';

=pod

=item new 

    my $writer = new Genome::Model::Tools::Pcap::Ace::Writer(\*STDOUT);    

=cut
sub new {
    croak("$pkg:new:no class given, quitting") if @_ < 1;
    my ($caller, $output) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self={};
    bless ($self, $class);
    
    %{$self->{'object_writers'}} = ( 
        assembly        => \&_write_assembly,
        read_position    => \&_write_read_position,
        contig          => \&_write_contig,
        'read'          => \&_write_read,
        assembly_tag     => \&_write_assembly_tag,
        contig_tag       => \&_write_contig_tag,
        read_tag         => \&_write_read_tag,
        base_segment     => \&_write_base_segment,
    );


    $self->{'output'} = $output;
    $self->{'width'} = 50;

    return $self;
}

=pod

=item set_output 

$writer->output(FH);

$FH = $writer->output();

=cut
sub output {
    my ($self,$FH) = @_;
    if (defined $FH) {
        $self->{'output'} = $FH;
    }
    return $self->{'output'};
}

sub width {
    my ($self,$width) = @_;
    if (defined $width) {
        $self->{'width'} = $width;
    }
    return $self->{'width'};
}

=pod

=item write_object 

    %contig = ( type => 'contig',
                name => 'Contig1',
                base_count => length $consensus,
                read_count => $read_count,
                base_seg_count => $base_seg_count,
                u_or_c => 'U',
                consensus => $consensus,
                base_qualities => \@base_qualities,
              );
    $writer->write_object(\%contig);

=cut
sub write_object {
    my ($self, $obj) = @_;
    my $OUT = $self->{'output'};
    my $type = $obj->{'type'};
    if (defined $type) {
        $self->{'object_writers'}->{$type}->($self,$OUT,$obj);
        if ($type eq 'base_segment') {
            $self->{'bs_last'} = 1;
        }
        else {
            $self->{'bs_last'} = 0;
        }
    }
    else {
        die "Undefined type in write_object\n";
    }
}

sub _write_assembly {
    my ($self, $OUT, $assembly) = @_;
    print $OUT "AS $assembly->{'contig_count'} $assembly->{'read_count'}\n\n"
}

sub _write_contig {
    my ($self, $OUT, $contig) = @_;
    print $OUT "\nCO $contig->{'name'} $contig->{'base_count'} $contig->{'read_count'} $contig->{'base_seg_count'} $contig->{'u_or_c'}\n";

    my $consensus = $contig->{'consensus'};
    $self->_write_sequence($OUT, $consensus);#, $contig->{'pads'});
    print $OUT "\n\nBQ";

    my $width = $self->width();
    my @bq = @{$contig->{'base_qualities'}};
    for (my $i = 0; $i < @bq; $i += 1) {
        if ($i % $width == 0) {
            print $OUT "\n";
        }
        print $OUT " $bq[$i]";
    }
    print $OUT "\n\n";
}

sub _write_sequence {
    my ($self, $OUT, $sequence, $pads) = @_;
    $sequence = $self->pad_sequence($sequence, $pads);
    my $seq_len = length($sequence);
    my $width = $self->width();
    for (my $i = 0;$i < $seq_len; $i += $width ) {
        print $OUT substr($sequence, $i, $i + ($width-1) < $seq_len ? $width : $seq_len - $i) . "\n";
    }
}

sub pad_sequence {
    my ($self, $sequence, $pads) = @_;
    my $pad_count = 0;
    foreach my $pad (@$pads) {
        my $pos = $pad->[0];
        my $pad_str = "*" x $pad->[1];
        $sequence = substr($sequence, 0, $pos + $pad_count) . $pad_str . substr($sequence, $pos + $pad_count);
        $pad_count += $pad->[1];
    }
    return $sequence;
}

sub _write_read_position {
    my ($self, $OUT, $read_pos) = @_;
    print $OUT "AF $read_pos->{'read_name'} $read_pos->{'u_or_c'} $read_pos->{'position'}\n";
}

sub _write_base_segment {
    my ($self, $OUT, $base_segment) = @_;
    print $OUT "BS $base_segment->{'start_pos'} $base_segment->{'end_pos'} $base_segment->{'read_name'}\n";
}

sub _write_read {
    my ($self, $OUT, $read) = @_;
    print $OUT "\nRD $read->{'name'} $read->{'padded_base_count'} $read->{'info_count'} $read->{'tag_count'}\n";
    my $sequence = $read->{'sequence'};
    $self->_write_sequence($OUT, $sequence);#, $read->{'pads'});
    print $OUT "\n";
    print $OUT "QA $read->{'qual_clip_start'} $read->{'qual_clip_end'} $read->{'align_clip_start'} $read->{'align_clip_end'}\n";
    print $OUT "DS";
    my %desc = %{$read->{'description'}};
    
    #CHROMAT_FILE PHD_FILE CHEM DYE TIME
    foreach my $key (qw(CHROMAT_FILE PHD_FILE CHEM DYE)) {
        if (exists $desc{$key}) {
            print $OUT " $key: $desc{$key}"; 
            delete $desc{$key};
        }
    }

    foreach my $key (sort keys %desc) {
        print $OUT " $key: $desc{$key}"; 
    }
    print $OUT "\n";
}

sub _write_assembly_tag {
    my ($self, $OUT, $tag) = @_;
    print $OUT "\nWA{\n$tag->{'tag_type'} $tag->{'program'} $tag->{'date'}\n$tag->{'data'}}\n";
}

sub _write_contig_tag {
    my ($self, $OUT, $tag) = @_;
    print $OUT "\nCT{\n$tag->{'contig_name'} $tag->{'tag_type'} $tag->{'program'} $tag->{'start_pos'} $tag->{'end_pos'} $tag->{'date'}";
    if ($tag->{'no_trans'}) {
        print $OUT " $tag->{no_trans}\n";
    }
    else {
        print $OUT "\n";
    }
    if (exists $tag->{'data'} && defined($tag->{'data'}) > 0) {
	#data already contain new line
        print $OUT "$tag->{'data'}";
    }
    print $OUT "}\n";
}

sub _write_read_tag {
    my ($self, $OUT, $tag) = @_;
    print $OUT "\nRT{\n$tag->{'read_name'} $tag->{'tag_type'} $tag->{'program'} $tag->{'start_pos'} $tag->{'end_pos'} $tag->{'date'}\n";
    print $OUT "$tag->{'data'}" if defined $tag->{data};
    print $OUT "}\n";
}

1;

#$Header$
