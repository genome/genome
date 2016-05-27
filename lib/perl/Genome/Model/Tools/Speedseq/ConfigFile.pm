package Genome::Model::Tools::Speedseq::ConfigFile;

use strict;
use warnings;

use Genome;
use File::Spec;

class Genome::Model::Tools::Speedseq::ConfigFile {
    is => 'Genome::SoftwareResult::StageableSimple',
    has => [
        'reference_sequence_build' => {
            is => 'Genome::Model::Build::ReferenceSequence',
            is_input => 1,
        },
        'speedseq_version' => {
            is => 'String',
            doc => 'The version of Speedseq used.',
            is_param => 1,
        },
    ],
};

sub _run {
    my $self = shift;

    my $per_chr = $self->reference_sequence_build->per_chromosome_fastas();
    my $fasta_dir = $per_chr->fasta_directory;
    my $old_config_path = Genome::Model::Tools::Speedseq::Base->config_path_for_version($self->speedseq_version);
    my $new_config_path = File::Spec->join($self->temp_staging_directory, 'speedseq.config');

    my $cmd = "sed 's|CNVNATOR_CHROMS_DIR=.*\$|CNVNATOR_CHROMS_DIR=$fasta_dir|' $old_config_path > $new_config_path";

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$old_config_path],
        output_files => [$new_config_path],
    );
    return $self;
}

sub config_file_path {
    my $self = shift;

    my ($allocation) = $self->disk_allocations;
    return File::Spec->join($allocation->absolute_path,'speedseq.config');
}


1;
