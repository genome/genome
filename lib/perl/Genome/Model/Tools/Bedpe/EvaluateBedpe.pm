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
            doc => 'The amount of slop (in b.p.). to be added to one set of breakpoints',
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
    $self->rawstats->{true_positive} = $self->_get_stat($self->bedpe, $self->gold_bedpe, 'both');
    $self->rawstats->{false_positive} = $self->_get_stat($self->bedpe, $self->gold_bedpe, 'notboth');
    $self->rawstats->{false_negative} = $self->_get_stat($self->gold_bedpe, $self->bedpe, 'notboth');
    $self->_set_derivative_stats;
    $self->print_stats;
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

sub _get_stat {
    my ($self, $file_a, $file_b, $type) = @_;

    my $output_file = Genome::Sys->create_temp_file_path;

    Genome::Model::Tools::BedTools::PairToPair->execute(
        $self->common_params,
        output_file => $output_file,
        input_file_a => $file_a,
        input_file_b => $file_b,
        intersection_type => $type,
    );
    return Genome::Sys->line_count($output_file);
}

sub _set_derivative_stats {
    my $self = shift;
    my $tp = $self->rawstats->{true_positive};
    my $fp = $self->rawstats->{false_positive};
    my $fn = $self->rawstats->{false_negative};

    $self->rawstats->{ppv} = $tp/($tp + $fp);
    $self->rawstats->{sensitivity} = $tp/($tp + $fn);
    $self->rawstats->{f1} = 2*$tp/(2*$tp + $fp + $fn);
}

sub print_stats {
    my $self = shift;

    for my $stat (sort keys %{$self->rawstats}) {
        printf("%s\t%s\n", $stat, $self->rawstats->{$stat});
    }
}
1;

