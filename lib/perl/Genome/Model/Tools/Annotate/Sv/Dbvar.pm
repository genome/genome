package Genome::Model::Tools::Annotate::Sv::Dbvar;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Annotate::Sv::Dbvar {
    is => "Genome::Model::Tools::Annotate::Sv::Base",
    has_input => [
        annotation_file => {
            is => 'Text',
            doc => 'File containing UCSC table',
            default => "/gsc/scripts/share/BreakAnnot_file/human_build37/GRCh37.remap.all.germline.ucsc.gff",
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

sub process_breakpoint_list {
    my ($self, $breakpoints_list) = @_;
    my %output;
    my $dbvar_annotation = $self->read_dbvar_annotation($self->annotation_file);
    $self->annotate_interval_matches($breakpoints_list, $dbvar_annotation, $self->breakpoint_wiggle_room, "dbvar_annotation", "bpB");
    foreach my $chr (keys %{$breakpoints_list}) {
        foreach my $item (@{$breakpoints_list->{$chr}}) {
            my $key = $self->get_key_from_item($item);
            $output{$key} = [$self->get_var_annotation($item, $item->{dbvar_annotation}->{bpB})];
        }
    }
    return \%output;
}

sub column_names {
    return ('dbvar');
}

sub read_dbvar_annotation {
    my $self = shift;
    my $file = shift;

    my %annotation;
    my $in = Genome::Sys->open_file_for_reading($file) || die "Unable to open annotation: $file\n";

    while (my $line = <$in>) {
       chomp $line;
       next if ($line =~ /^\#/);
       my $p;
       my @fields = split /\t+/, $line;
       $p->{chrom} = $fields[0];
       $p->{chromStart} = $fields[3];
       $p->{chromEnd} = $fields[4];
       $p->{name} = $fields[8];
       $p->{chrom} =~ s/chr//;
       push @{$annotation{$p->{chrom}}{$p->{chromEnd}}{$p->{chromStart}}}, $p;
    }
    return \%annotation;
}

1;

