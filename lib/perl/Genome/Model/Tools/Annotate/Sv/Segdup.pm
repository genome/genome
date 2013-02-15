package Genome::Model::Tools::Annotate::Sv::Segdup;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Annotate::Sv::Segdup {
    is => "Genome::Model::Tools::Annotate::Sv::IntervalAnnotator",
    has_input => [
        annotation_file => {
            is => 'String',
            doc => 'File containing UCSC table',
            default => "/gsc/scripts/share/BreakAnnot_file/human_build37/Human.Feb2009.SegDups.tab",
        }
    ],
};

sub process_breakpoint_list {
    my ($self, $breakpoints_list) = @_;
    my %output;
    my $segdup_annotation = $self->read_ucsc_annotation($self->annotation_file);
    $self->annotate_interval_matches($breakpoints_list, $segdup_annotation, $self->breakpoint_wiggle_room, "segdup_annotation");
    foreach my $chr (keys %{$breakpoints_list}) {
        foreach my $item (@{$breakpoints_list->{$chr}}) {
            my $key = $self->get_key_from_item($item);
            $output{$key} = $item->{segdup_annotation};
            if (defined $item->{segdup_annotation}) {
                my @segdup = map {$_->{name}} @{$item->{segdup_annotation}};
                $output{$key} = [join(",", @segdup)];
            }
            else {
                $output{$key} = ["N/A"];
            }
        }
    }
    return \%output;
}

sub column_names {
    return ('segdup');
}

1;

