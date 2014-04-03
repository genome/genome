package Genome::Model::Tools::Cufflinks::Assemble;

use strict;
use warnings;

use Genome;
use Cwd;
use version;

class Genome::Model::Tools::Cufflinks::Assemble {
    is => 'Genome::Model::Tools::Cufflinks',
    has_input => [
        input_file => {
            doc => 'The sam/bam file to generate transcripts and expression levels from.',
        },
        params => {
            doc => 'Any additional parameters to pass to cufflinks',
            is_optional => 1,
        },
        output_directory => { doc => 'The directory to write all output files to', is => 'Text',},
    ],
    has_output => [
        transcripts_file => {
            is_optional => 1,
        },
        transcript_expression_file => {
            is_optional => 1,
        },
        gene_expression_file => {
            is_optional => 1,
        },
        assembler_output_file => {
            is_optional => 1,
        },
        gene_fpkm_file => {
            is_optional => 1,
        },
        isoform_fpkm_file => {
            is_optional => 1,
        },
    ],
};

sub command {
    my $self = shift;
    my ($params) = @_;

    my $cmd = $self->cufflinks_path
        . ' --no-update-check '
        . $params .' '
        . $self->input_file
        . ' > '. $self->assembler_output_file
        . ' 2>&1';

    return $cmd;
}

sub execute {
    my $self = shift;

    my $output_directory = $self->output_directory;

    my $params = $self->params || '';

    $self->transcripts_file($output_directory .'/transcripts.gtf');
    $self->assembler_output_file($output_directory .'/cufflinks.out');
    my @output_files = ($self->transcripts_file, $self->assembler_output_file);

    my $cwd;
    # Versions prior to 1.0.0 wrote some output files to cwd
    if ( version->parse($self->use_version) < version->parse('1.0.0') ) {
        $cwd = getcwd;
        unless (chdir($output_directory)) {
            $self->error_message('Failed to change cwd to '. $output_directory);
            die($self->error_message);
        }

        $self->transcript_expression_file($output_directory .'/transcripts.expr');
        $self->gene_expression_file($output_directory .'/genes.expr');

        push @output_files, $self->transcript_expression_file;
        push @output_files, $self->gene_expression_file;
    } else {
        $params .= ' -o '. $output_directory;

        $self->gene_fpkm_file($output_directory .'/genes.fpkm_tracking');
        $self->isoform_fpkm_file($output_directory .'/isoforms.fpkm_tracking');

        push @output_files, $self->gene_fpkm_file;
        push @output_files, $self->isoform_fpkm_file;
    }

    if (version->parse($self->use_version) >= version->parse('0.9.0')) {
        # The progress bar since v0.9.0 is causing massive(50MB) log files
        # TODO: should we parse the params to see if -q or -v are already defined?
        $params .= ' -q ';
    }

    Genome::Sys->shellcmd(
        cmd => $self->command($params),
        input_files => [$self->input_file],
        output_files => \@output_files,
    );
    if (defined($cwd)) {
        unless (chdir($cwd)) {
            $self->error_message('Failed to cd back to original directory: '. $cwd);
            die($self->error_message);
        }
    }
    return 1;
}


1;
