package Genome::Model::PhenotypeCorrelation::Command::ClinicalCorrelationWrapper;

use File::Basename qw/basename/;
use Data::Dumper;
use Genome;
use Carp qw/confess/;
use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::ClinicalCorrelationWrapper {
    is => ["Genome::Command::Base"],
    doc => "Run the clinical correlation tool from the music suite.",
    has => [
        variant_matrix => {
            is => "File",
            doc => "Variant matrix to analyze",
            is_input => 1,
        },
        output_directory => {
            is => "File",
            doc => "Where to store the output files",
            is_input => 1,
        },
        sample_list_file => {
            is => "String",
            doc => 'File containing samples names, 1 per line, for input into MuSiC',
            is_input => 1,
        },
        clinical_data_file => {
            is => "String",
            doc => 'File containing clinical data',
            is_input => 1,
        },
        categorical_clinical_data_file => {
            is => "String",
            doc => "File containing categorical clinical data for Fisher's exact test",
            is_input => 1,
            is_optional => 1,
        },
        glm_model_file => {
            is => "String",
            doc => 'File containing the model specification for this analysis',
            is_input => 1,
        },
    ],
};

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
This command exists in order to facilitate parallelization of
music's clinical correlation command. The output file will
be named \$output_directory/\$variant_matrix.results.glm.tsv.
This allows us to split a variant matrix into many parts and
process them in parallel easily with Workflow.
EOS
}

sub execute {
    my $self = shift;

    my $out_dir = $self->output_directory;
    my $variant_matrix = $self->variant_matrix;
    my $glm_model_file = $self->glm_model_file;
    my $clin_file = $self->clinical_data_file;
    my $categorical_file = $self->categorical_clinical_data_file;
    my $samples_file = $self->sample_list_file;
    # the .glm.tsv is added by the actual ClinicalCorrelation command
    my $output_file = join("/", $out_dir, basename($variant_matrix)) . ".results";

    my %params = (
        genetic_data_type => "variant",
        bam_list => $samples_file,
        output_file => $output_file,
        glm_model_file => $glm_model_file,
        glm_clinical_data_file => $clin_file,
        input_clinical_correlation_matrix_file => $variant_matrix,
    );
    $params{categorical_clinical_data_file} = $categorical_file if $categorical_file;

    $self->status_message("Preparing to run clinical correlation with params:");
    $self->status_message(Dumper(\%params));

    my $command = Genome::Model::Tools::Music::ClinicalCorrelation->create(
        %params
    );
    confess "Clinical correlation failed for $variant_matrix!" unless $command->execute;
    $self->status_message("Output files:\n\t" . join("\n\t", glob("$output_file*")));

    return 1;
};

1;
