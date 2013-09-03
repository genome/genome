package Genome::Model::Tools::Dindel::AnalyzeWindowFile;

use strict;
use warnings;

use Genome;

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
        input_bam => {
            is => 'Path',
        },
        ref_fasta => {
            is => 'Path',
        },
    ],
    has_calculated_output => [
        output_bam => {
            is => 'Path',
            calculate =>  q{ $output_prefix . "_realigned.merged.bam"},
            calculate_from => ['output_prefix'],
        },
        output_log => {
            is => 'Path',
            calculate =>  q{ $output_prefix . "_realigned.windows.log"},
            calculate_from => ['output_prefix'],
        },
        output_glf => {
            is => 'Path',
            calculate =>  q{ $output_prefix . ".glf.txt"},
            calculate_from => ['output_prefix'],
        },
    ],
    has_optional_transient => {
        output_prefix => {
            is_calculated => 1,
            calculate =>  q{ File::Spec->join($output_directory, "dindel") },
            calculate_from => ['output_directory'],
        },
    },
};

sub help_brief {
    'Actual slow part of dindel-- analysis'
}

sub execute {
    my $self = shift;
    $self->create_output_directory();

    my @cmd = (
        $self->dindel_executable,
        '--analysis', 'indels',
        '--doDiploid',
        '--bamFile', $self->input_bam,
        '--varFile', $self->window_file,
        '--outputFile', $self->output_prefix,
        '--ref', $self->ref_fasta,
        '--libFile', $self->library_metrics_file,
        '--outputRealignedBAM',
        '--processRealignedBAM', $self->merge_bam_script,
    );
    my @input_files = (
        $self->input_bam,
        $self->window_file,
        $self->ref_fasta,
        $self->merge_bam_script,
    );
    my $rv = $self->shellcmd_arrayref(
        cmd => \@cmd,
        input_files => \@input_files,
    );
    unless ($rv) {
        die "Failed to complete dindel analysis.";
    }

    $self->convert_sam_to_bam();
    return 1;
}

sub merge_bam_script {
    my $self = shift;

    my (undef, $directory, undef) = File::Spec->splitpath(__FILE__);
    return File::Spec->join($directory, 'merge_bam_callback.pl');
}

sub convert_sam_to_bam {
    my $self = shift;

    my $reheadered_sam_file = $self->reheader_sam_file();
    Genome::Model::Tools::Sam::SamToBam->execute(
        sam_file => $reheadered_sam_file,
        keep_sam => 0,
        bam_file => $self->output_bam,
        is_sorted => 0,
        index_bam => 0,
    );
}

sub reheader_sam_file {
    my $self = shift;

    (my $output_sam_file = $self->output_bam) =~ s/bam$/sam/;

    my $reheadered_sam_file = $self->output_prefix . ".reheadered.sam";
    my $cmd = "cat <(samtools view -H %s) %s > %s ";
    my @input_files = (
        $self->input_bam,
        $output_sam_file,
    );
    my @output_files = (
        $reheadered_sam_file,
    );
    my $rv = Genome::Sys->shellcmd(
        cmd => sprintf($cmd, @input_files, @output_files),
        input_files => \@input_files,
        output_files => \@output_files,
    );
    unless ($rv) {
        die "Failed to reheader sam file $output_sam_file";
    }

    unlink $output_sam_file;
    return $reheadered_sam_file
}

1;
