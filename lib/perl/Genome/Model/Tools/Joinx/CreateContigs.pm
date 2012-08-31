package Genome::Model::Tools::Joinx::CreateContigs;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Joinx::CreateContigs {
    is => 'Genome::Model::Tools::Joinx',
    has_input => [
        reference => {
            is => 'Text',
            doc => 'Reference sequence fasta file',
            shell_args_position => 1,
        },
        variants => {
            is => 'Text',
            doc => 'Variants bed file',
            shell_args_position => 2,
        },
        output_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'output fasta file of contigs (default=stdout)',
        },
        flank => {
            is => 'Number',
            doc => 'Amount of flanking to pad contigs with on either side',
            default_value => 99,
        },
        quality => {
            is => 'Number',
            doc => 'Minimum quality of variants to consider',
            default_value => 0,
        },
    ],
};

sub help_brief {
    "Generate contigs given a reference sequence and bed file of variants.";
}

sub help_synopsis {
    my $self = shift;
    "gmt joinx create-contigs ref.fa variants.bed --flank 75 --quality 5";
}

sub execute {
    my $self = shift;
    my $output = $self->output_file || "-";
    my @args = (
        $self->joinx_path,
        "create-contigs",
        "--reference", $self->reference,
        "--variants", $self->variants,
        "--output-file", $output,
        "--flank", $self->flank,
        "--quality", $self->quality,
    );
    my $cmd = join(" ", @args);

    my %params = (
        cmd => $cmd,
        allow_zero_size_output_files=>1,
    );
    $params{output_files} = [$output] if $output ne "-";
    Genome::Sys->shellcmd(%params);

    return 1;
}

1;
