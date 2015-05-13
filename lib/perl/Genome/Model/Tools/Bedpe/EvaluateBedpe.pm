package Genome::Model::Tools::Bedpe::EvaluateBedpe;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Bedpe::EvaluateBedpe {
    is => 'Command::V2',
    has_input => [
        bedpe => {
            is => 'File',
            doc => 'bedpe file to be evaluated',
        },
        gold_bedpe => {
            is => 'Path',
            doc => 'bedpe file with gold standard breakpoints'
        },
        bedtools_version => {
            is => 'Text',
        },
        slop => {
            is => 'Integer',
            doc => 'The amount of slop (in b.p.). to be added to each footprint of A. *Note*: Slop is subtracted from start1 and start2 and added to end1 and end2.',
        },
    ],
    has_transient_optional_output => [
        rawstats => {
            is => "HASH",
            doc => "The raw stats generated during primary execution",
        },
    ],
};

sub execute {
    my $self = shift;
    $self->rawstats({});

    my $true_positives_file = Genome::Sys->create_temp_file_path;
    Genome::Model::Tools::BedTools::PairToPair->execute(
        $self->common_params,
        output_file => $true_positives_file,
        input_file_a => $self->bedpe,
        input_file_b => $self->gold_bedpe,
        intersection_type => 'both',
    );
    $self->rawstats->{true_positive} = Genome::Sys->line_count($true_positives_file);

    my $false_positives_file = Genome::Sys->create_temp_file_path;
    Genome::Model::Tools::BedTools::PairToPair->execute(
        $self->common_params,
        output_file => $false_positives_file,
        input_file_a => $self->bedpe,
        input_file_b => $self->gold_bedpe,
        intersection_type => 'notboth',
    );
    $self->rawstats->{false_positive} = Genome::Sys->line_count($false_positives_file);

    my $false_negatives_file = Genome::Sys->create_temp_file_path;
    Genome::Model::Tools::BedTools::PairToPair->execute(
        $self->common_params,
        output_file => $false_negatives_file,
        input_file_a => $self->gold_bedpe,
        input_file_b => $self->bedpe,
        intersection_type => 'notboth',
    );
    $self->rawstats->{false_negative} = Genome::Sys->line_count($false_negatives_file);
    return 1;
}

sub common_params {
    my $self = shift;
    return (
        slop => $self->slop,
        slop_strand => "+-",
        ignore_strand => 1,
        require_different_names => 0,
        use_version => $self->bedtools_version,
    );
}

1;

