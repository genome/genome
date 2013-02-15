package Genome::Model::Tools::Annotate::Sv::Base;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Annotate::Sv::Base{
    is => "Command::V2",
};

sub process_breakpoint_list {
    #override in subclass;
    #interface for sv annotators
};

sub column_names {
    #override in subclass;
    #interface for sv annotator
};

sub get_key_from_item {
    my $class = shift;
    my $item = shift;
    return join("--", $item->{chrA}, $item->{bpA}, $item->{chrB}, $item->{bpB}, $item->{event});
}

sub annotate_interval_matches {
    #both breakpoints need to match within some wiggle room
    my $self = shift;
    my $positions = shift;
    my $annotation = shift;
    my $annot_length = shift;
    my $tag = shift;

    foreach my $chr (keys %$positions) {
        my @sorted_items = sort {$a->{bpB}<=>$b->{bpB}} (@{$positions->{$chr}});
        my @sorted_positions = map{$_->{bpB}} @sorted_items;
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
                            if ($position_item->{bpB} eq $pos) {
                                push @{$position_item->{$tag}}, $var;
                            }
                        }
                    }
                }
            }
        }
    }
    return 1;
}

sub read_ucsc_annotation{
    my ($self, $file) = @_;
    my %annotation;
    open (ANNOTATION, "<$file") || die "Unable to open annotation: $file\n";
    while (<ANNOTATION>) {
        chomp;
        next if /^\#/;
        my $p;
        my @extra;
        ($p->{bin},$p->{chrom},$p->{chromStart},$p->{chromEnd},$p->{name},@extra) = split /\t+/;
        $p->{chrom} =~ s/chr//;
        $p->{extra} = \@extra;
        push @{$annotation{$p->{chrom}}{$p->{chromEnd}}{$p->{chromStart}}}, $p;
    }
    close ANNOTATION;
    return \%annotation;
}

1;
