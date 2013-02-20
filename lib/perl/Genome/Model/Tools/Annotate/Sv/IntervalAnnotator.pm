package Genome::Model::Tools::Annotate::Sv::IntervalAnnotator;

use strict;
use warnings;
use Genome;
use List::Util qw(min max);

class Genome::Model::Tools::Annotate::Sv::IntervalAnnotator {
    is => "Genome::Model::Tools::Annotate::Sv::Base",
    has_input => [
        annotation_file => {
            is => 'Text',
            doc => 'File containing UCSC table',
        },
        breakpoint_wiggle_room => {
            is => 'Number',
            doc => 'Distance between breakpoint and annotated breakpoint within which they are considered the same, in bp',
            default => 200,
        },
        overlap_fraction => {
            is => 'Number',
            doc => 'Fraction of overlap (reciprocal) required to hit',
            default => 0.5,
        },
    ],
};

#TODO NEEDS to be rewritten - I don't think it is right.
#First of all, it stops at one end position rather than getting all intervals that cross the breakpoint.
#Second of all, it only considers the 2nd breakpoint of the SV, not the first breakpoint.  I think
#annotations crossing either breakpoint should probably be considered.
sub annotate_interval_matches {
    #both breakpoints need to match within some wiggle room
    my $self = shift;
    my $positions = shift;
    my $annotation = shift;
    my $annot_length = shift;
    my $tag = shift;
    my $breakpoint_key = shift;

    foreach my $chr (keys %$positions) {
        my @sorted_items = sort {$a->{$breakpoint_key}<=>$b->{$breakpoint_key}} (@{$positions->{$chr}});
        my @sorted_positions = map{$_->{$breakpoint_key}} @sorted_items;
        my %annotated_output;

        my @chromEnds = sort {$a<=>$b} keys %{$annotation->{$chr}};
        for my $pos (@sorted_positions) {
            while (@chromEnds>0 && $pos>$chromEnds[0]+$annot_length) {
                shift @chromEnds;
            }
            next unless @chromEnds>0;
            for my $start (keys %{$$annotation{$chr}{$chromEnds[0]}}) {
                if ($pos>=$start-$annot_length) {
                    for my $var (@{$$annotation{$chr}{$chromEnds[0]}{$start}}){
                        foreach my $position_item (@{$positions->{$chr}}) {
                            if ($position_item->{$breakpoint_key} eq $pos) {
                                push @{$position_item->{$tag}->{$breakpoint_key}}, $var;
                            }
                        }
                    }
                }
            }
        }
    }
    return 1;
}

sub get_var_annotation {
    my ($self, $item, $annotation_ref) = @_;
    
    my $varreport = "N/A";
    my @vars;
    my $frac = $self->overlap_fraction;
    
    if (defined $annotation_ref) {
        foreach my $var (@$annotation_ref) {
            my $pos1 = min($item->{bpB}, $var->{chromEnd});
            my $pos2 = max($item->{bpA}, $var->{chromStart});
            my $overlap = $pos1-$pos2+1;
            my $ratio1 = $overlap/(abs($item->{bpB}-$item->{bpA})+1);
            my $ratio2 = $overlap/(abs($var->{chromEnd}-$var->{chromStart})+1);
            if ($ratio1 >= $frac && $ratio2 >= $frac ) {
                push @vars, $var;
            }
        }

        if (@vars) {
            $varreport = join(",",map{$_->{name}} @vars);
        }
    }
    return $varreport;
}

sub read_ucsc_annotation{
    my ($self, $file) = @_;
    my %annotation;
    my $in = Genome::Sys->open_file_for_reading($file) || die "Unable to open annotation: $file\n";
    while (<$in>) {
        chomp;
        next if /^\#/;
        my $p;
        my @extra;
        ($p->{bin},$p->{chrom},$p->{chromStart},$p->{chromEnd},$p->{name},@extra) = split /\t+/;
        $p->{chrom} =~ s/chr//;
        $p->{extra} = \@extra;
        push @{$annotation{$p->{chrom}}{$p->{chromEnd}}{$p->{chromStart}}}, $p;
    }
    $in->close;
    return \%annotation;
}

1;

