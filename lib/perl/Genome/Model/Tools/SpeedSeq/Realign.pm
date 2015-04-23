package Genome::Model::Tools::SpeedSeq::Realign;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SpeedSeq::Realign {
    is => 'Genome::Model::Tools::SpeedSeq::AlignBase',
    has_input => [
        bams => {
            is => 'Text',
            doc => 'BAM file(s) (must contain read group tags)',
            is_many => 1,
            example_values => ['in.bam'],
        },
    ],
    has_param => [
        rename_reads => {
            is => 'Boolean',
            doc => 'rename reads for smaller file size',
            is_optional => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    # OPTIONS aka params
    my $options = $self->_resolve_options_string();
    if ($self->rename_reads) {
        $options .= ' -n';
    }

    # INPUTS
    my @input_files;
    my $inputs_string = $self->reference_fasta;
    push @input_files, $self->reference_fasta;
    for my $bam ($self->bams) {
        push @input_files, $bam;
        $inputs_string .= ' '. $bam;
    }
    if ($self->config_file) {
        push @input_files, $self->config_file;
    }

    # COMMAND
    my $cmd = $self->speedseq_path .' realign '. $options .' '. $inputs_string;

    # OUTPUTS
    my $basename = 'in.realign';
    if ($self->output_prefix) {
        $basename = $self->output_prefix;
    }
    my @output_files = $self->_resolve_output_file_paths($basename);
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => \@input_files,
        output_files => \@output_files,
    );
    $self->output_files(\@output_files);

    return 1;
}

1;
