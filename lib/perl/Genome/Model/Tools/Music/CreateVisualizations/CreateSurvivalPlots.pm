package Genome::Model::Tools::Music::CreateVisualizations::CreateSurvivalPlots;

use warnings;
use strict;
use Genome;

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::CreateVisualizations::CreateSurvivalPlots {
    is => 'Command::V2',
    has_input => [
        output_dir => {is => 'Text', doc => 'Output directory path'},
        maf_file => {is => 'Text', doc => 'final maf output file'},
        categorical_clinical_data_file => {is => 'Text', is_optional => 1, doc => 'Table of samples (y) vs. categorical clinical data category (x)'},
        numeric_clinical_data_file => {is => 'Text', is_optional => 1, doc => 'Table of samples (y) vs. numeric clinical data category (x)'},
        bam_list => {is => 'Text', doc => 'Tab delimited list of BAM files [sample_name normal_bam tumor_bam]'},
    ],
    
};

sub execute {
    my $self = shift;
    my $cmd = Genome::Model::Tools::Music::Survival->create(
        output_dir => $self->output_dir,
        maf_file => $self->maf_file,
        bam_list => $self->bam_list,
        categorical_clinical_data_file => $self->categorical_clinical_data_file,
        numeric_clinical_data_file => $self->numeric_clinical_data_file,
    );
    my $rv = eval{$cmd->execute()};
    if($@){
        my $error = $@;
        $self->error_message('Error running ' . $cmd->command_name . ': ' . $error);
        return;
    }
    return 1;
}

1;
