package Genome::Model::Tools::Annotate::Sv::Dbsnp;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Annotate::Sv::Dbsnp {
    is => "Genome::Model::Tools::Annotate::Sv::IntervalAnnotator",
    has_input => [
        annotation_file => {
            is => 'Text',
            doc => 'File containing UCSC table',
            default => "/gsc/scripts/share/BreakAnnot_file/human_build37/dbsnp132.indel.named.csv",
        },
    ],
};

sub process_breakpoint_list {
    my ($self, $breakpoints_list) = @_;
    my %output;
    my $dbsnp_annotation = $self->read_ucsc_annotation($self->annotation_file);
    $self->annotate_interval_matches($breakpoints_list, $dbsnp_annotation, $self->breakpoint_wiggle_room, "dbsnp_annotation");
    foreach my $chr (keys %{$breakpoints_list}) {
        foreach my $item (@{$breakpoints_list->{$chr}}) {
            my $key = $self->get_key_from_item($item);
            $output{$key} = $item->{dbsnp_annotation};
            if (defined $item->{dbsnp_annotation}) {
                my @dbsnp = map {$_->{name}} @{$item->{dbsnp_annotation}};
                $output{$key} = [join(",", @dbsnp)];
            }
            else {
                $output{$key} = ["N/A"];
            }
        }
    }
    return \%output;
}

sub column_names {
    return ('dbsnp');
}

1;

