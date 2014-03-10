package Genome::Model::PhenotypeCorrelation::Command::BurdenTestWrapper;

use Genome;
use Data::Dumper;
use JSON;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::BurdenTestWrapper {
    is => "Command::V2",
    has_input => [
        params => {
            is => "Text",
            doc => "JSON encoded hashref containing gene, phenotype, analysis_data_type, covariates",
        },
        option_file => {
            is => 'Text',
            doc => 'File specifying informatin for the burden test to run',
        },
        permutations => {
            is => 'Integer',
            doc => 'The number of permutations to perform for calculating the p-value',
            default => 10000,
        },
        seed => {
            is => 'Integer',
            doc => 'The seed for the random number generator used to do the permutations',
            default => 123,
        },
        p_value_permutations => {
            is => 'Integer',
            doc => 'The number of permutations to perform on the p-values',
            default => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    my $json = new JSON;
    my $params = $json->decode($self->params);
    $self->status_message("Decoded params: " . Dumper($params));
    my %cmd_params = (
        analysis_data_type => $params->{analysis_data_type},
        phenotype_name => $params->{phenotype},
        gene_name => $params->{gene},
        covariates => join("+", @{$params->{covariates}}),
        option_file => $self->option_file,
        permutations => $self->permutations,
        seed => $self->seed,
        p_value_permutations => $self->p_value_permutations,
    );

    printf "Params: %s\n", Dumper(\%cmd_params);

    my $cmd = Genome::Model::Tools::Germline::BurdenTest->create(%cmd_params);
    return $cmd->execute;
}

1;
