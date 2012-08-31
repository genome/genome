package Genome::Model::Tools::DetectVariants::Somatic::Pindel::Large;

use warnings;
use strict;

use Genome;
use Workflow;
use File::Copy;
use Workflow::Simple;
use Cwd;

class Genome::Model::Tools::DetectVariants::Somatic::Pindel::Large {
    is => ['Genome::Model::Tools::DetectVariants::Somatic::Pindel'],
    doc => "Runs the pindel pipeline , provides 32GB lsf_resource",
    has_param => [
        lsf_resource => {
            default_value => "-M 32000000 -R 'select[type==LINUX64 && mem>32000] rusage[mem=32000]'",
        },
    ],
};

sub _run_converter {
    my $self = shift;
    my $throwaway = shift;
    my $converter = 'Genome::Model::Tools::Bed::Convert::Indel::PindelToBed';
    my $source = shift;
    
    my $output = $source . '.bed'; #shift; #TODO Possibly create accessors for the bed files instead of hard-coding this
    
    my $command = $converter->create(
        source => $source,
        output => $output, 
        reference_sequence_input => $self->reference_sequence_input,
    );
    
    unless($command->execute) {
        $self->error_message('Failed to convert ' . $source . ' to the standard format.');
        return;
    }

    return 1;
}

