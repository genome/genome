package Genome::Model::Tools::Annotate::Sv::Dbsnp;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Annotate::Sv::Dbsnp {
    is => "Genome::Model::Tools::Annotate::Sv::Base",
    doc => "Annotate overlaps with dbsnp SVs",
    has_input => [
        annotation_file => {
            is => 'Text',
            doc => 'File containing UCSC table',
            example_values => ["/gsc/scripts/share/BreakAnnot_file/human_build37/dbsnp132.indel.named.csv"],
        },
        overlap_fraction => {
            is => 'Number',
            doc => 'Fraction of overlap (reciprocal) required to hit',
            default => 0.5,
        },

    ],
};

sub help_detail {
    return "Determines whether the SV breakpoints match a dbSNP SV within some distance.  It also checks to see that the SV and the dbSNP SV reciprocally overlap each other by a given fraction.";
}

sub process_breakpoint_list {
    my ($self, $breakpoints_list) = @_;
    my %output;
    my $dbsnp_annotation = $self->read_ucsc_annotation($self->annotation_file);
    $self->annotate_interval_overlaps($breakpoints_list, $dbsnp_annotation, "dbsnp_annotation");
    foreach my $chr (keys %{$breakpoints_list}) {
        foreach my $item (@{$breakpoints_list->{$chr}}) {
            my $key = $self->get_key_from_item($item);
            my $a = $item->{dbsnp_annotation}->{bpA};
            my $b = $item->{dbsnp_annotation}->{bpB};
            
            $output{$key} = [$self->get_var_annotation($item, $a, $b)];
        }
    }
    return \%output;
}

sub column_names {
    return ('dbsnp');
}

1;

