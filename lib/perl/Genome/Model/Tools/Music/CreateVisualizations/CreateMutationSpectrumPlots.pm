package Genome::Model::Tools::Music::CreateVisualizations::CreateMutationSpectrumPlots;

use warnings;
use strict;
use Genome;

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::CreateVisualizations::CreateMutationSpectrumPlots {
    is => 'Command::V2',
    has_input => [
        output_dir => {is => 'Text', doc => 'Output directory path'},
        bed_list => {is => 'Text', doc => 'Tab delimited list of BED files [label bed_file]'},
    ],
};

sub execute {
    my $self = shift;
    my $output_file = join('/', $self->output_dir, 'mutation_spectrum_plot.pdf');
    my $mut_spec_file = join('/', $self->output_dir, 'mutation_spectrum_file.tsv');
    my $cmd = Genome::Model::Tools::Analysis::SummarizeMutationSpectrum->create(
        input_file => $self->bed_list, 
        output_file => $output_file,
        mut_spec_file => $mut_spec_file,
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
