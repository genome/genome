package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BackfillIndelVcfHelper;

use strict;
use warnings;

use above 'Genome';
use File::Spec;

class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BackfillIndelVcfHelper {
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
    has_calculated_output => [
        varscan_output_file => {
            is => 'Path',
            calculate =>  q{ File::Spec->join($output_directory, "varscan_consensus.vcf.gz") },
            calculate_from => ['output_directory'],
        },
        joinx_output_file => {
            is => 'Path',
            calculate =>  q{ File::Spec->join($output_directory, "joinx_normalized.vcf.gz") },
            calculate_from => ['output_directory'],
        },
        output_vcf => {
            is => 'Boolean',
            calculate =>  q{ 1; },
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
