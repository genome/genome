package Genome::Model::Tools::Dindel::AnalyzeWindowFile;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::Dindel::AnalyzeWindowFile {
    is => 'Genome::Model::Tools::Dindel::Base',
    has_input => [
        window_file => {
            is => 'Path',
        },
        library_metrics_file => {
            is => 'Path',
            doc => 'from step one... getCigarIndels',
        },
        bam_file => {
            is => 'Path',
        },
        output_prefix => {
            is => 'Path',
            is_output => 1,
        },
        ref_fasta => {
            is => 'Path',
        },
        output_bam => {
            is => 'Boolean',
            default => 0,
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Actual slow part of dindel-- analysis'
}

sub help_synopsis {
    return <<EOS
EOS
}

sub help_detail {
    return <<EOS
EOS
}


sub execute {
    my $self = shift;

    my @cmd = (
        $self->dindel_executable,
        '--analysis', 'indels',
        '--doDiploid',
        '--bamFile', $self->bam_file,
        '--varFile', $self->window_file,
        '--outputFile', $self->output_prefix,
        '--ref', $self->ref_fasta,
        '--libFile', $self->library_metrics_file,
    );

    if ($self->output_bam) {
        my (undef, $callback_script_path) = File::Basename::fileparse(__FILE__);
        my $callback_script = $callback_script_path . "merge_bam_callback.pl";
        unless(-s $callback_script) {
            $self->error_message("unable to locate dindel helper script at location: $callback_script\n");
            return undef;
        }

        push @cmd, '--outputRealignedBAM';
        push @cmd, '--processRealignedBAM', $callback_script;
    }


    return Genome::Sys->shellcmd_arrayref(
        cmd => \@cmd,
        input_files => [
            $self->bam_file,
            $self->window_file,
            $self->ref_fasta,
        ],
    );
}

1;
