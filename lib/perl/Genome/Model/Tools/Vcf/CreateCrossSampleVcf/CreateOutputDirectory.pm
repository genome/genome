package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateOutputDirectory;

use strict;
use warnings;

use above 'Genome';

class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateOutputDirectory {
    is => 'Command::V2',
    has_input => [
        sample_name => {
            is => 'Text',
        },
        base_directory => {
            is => 'Path',
        },
    ],
    has_optional_output => [
        output_directory => {
            is => 'Path',
        },
    ],
};

sub execute {
    my $self = shift;

    my $sample = Genome::Sample->get(name => $self->sample_name);
    unless ($sample) {
        die sprintf("Couldn't find sample with name: %s", $self->sample_name);
    }

    my $output_directory = File::Spec->join($self->base_directory, $sample->id);
    if (!-d $output_directory) {
        Genome::Sys->create_directory($output_directory);
    }
    unless (-d $output_directory) {
        die sprintf("Couldn't create directory: %s", $output_directory);
    }
    $self->output_directory($output_directory);

    return 1;
}

1;
