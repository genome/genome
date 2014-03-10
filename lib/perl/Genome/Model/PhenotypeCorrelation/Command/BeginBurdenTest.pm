package Genome::Model::PhenotypeCorrelation::Command::BeginBurdenTest;

use Genome;
use Genome::File::Vep::Reader;
use Data::Dumper;
use JSON;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::BeginBurdenTest {
    is => "Command::V2",
    has_input => [
        option_file => {
            is => 'Text',
            doc => "R option file",
        },
        genes => {
            is_many => 1,
            is => "Text",
            doc => "The list of gene names to process",
        },
        glm_model_file => {
            is => 'Text',
            doc => 'File outlining the type of model, response variable, covariants, etc. for the GLM analysis. (See DESCRIPTION).',
        },
    ],
    has_output => [
        job_params => {
            is => "Text",
            is_many => 1,
            doc => "list of JSON encoded hashrefs containing gene, phenotype, analysis_data_type, covariates",
        },
    ]
};

sub execute {
    my $self = shift;

    my $glm_config = Genome::Model::PhenotypeCorrelation::GlmConfig->from_file($self->glm_model_file);

    my @job_params;
    my $json = new JSON;
    for my $attr ($glm_config->attributes) {
        my $phenotype = $attr->{attr_name};
        my $analysis_data_type = $attr->{type};
        my $covariates = $attr->{covariates};
        for my $gene ($self->genes) {
            my $params = $json->encode({
                    gene => $gene,
                    phenotype => $phenotype,
                    analysis_data_type => $analysis_data_type,
                    covariates => $covariates,
                    });

            push(@job_params, $params);
        }
    }
    $self->job_params(\@job_params);
    return 1;
}

1;
