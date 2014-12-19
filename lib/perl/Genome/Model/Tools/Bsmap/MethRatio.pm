package Genome::Model::Tools::Bsmap::MethRatio;

use strict;
use warnings;
use Genome;

use Genome::Utility::Text qw(sanitize_string_for_filesystem);
use File::Spec;

my $DEFAULT_VERSION = '2.74';
my $METHRATIO_COMMAND = 'methratio.py';

my %METHRATIO_VERSIONS = (
    2.6 => File::Spec->join('/gscuser/cmiller/usr/src/bsmap-2.6', $METHRATIO_COMMAND),
    2.74 => File::Spec->join('/gsc/pkg/bio/bsmap/bsmap-2.74', $METHRATIO_COMMAND),
);

class Genome::Model::Tools::Bsmap::MethRatio {
    is => 'Command',
    has => [
        bam_file => {
            is => 'Text',
            doc => 'The bam file to do the counting on (must be a product of bsmap alignment)',
            is_input => 1,
        },
        output_file => {
            is => 'Text',
            doc => 'File name for methyl counts',
            is_output => 1,
            is_input => 1,
            default_value => 'snvs.hq',
        },
        output_directory => {
            is => 'Text',
            doc => 'Where to output the methyl counts',
            is_output => 1,
            is_input => 1,
        },
        reference => {
            is => 'Text',
            doc => '36, 37, or a path to the reference fasta',
            is_input => 1,
        },
        version => {
            is => 'Version',
            is_optional => 1,
            is_input => 1,
            default_value => $DEFAULT_VERSION,
            doc => 'Version of methratio to use',
        },
    ],
    has_optional => [
        chromosome => {
            is => 'Text',
            doc => 'process only this chromosome',
            is_input => 1,
        },
        output_zeros => {
            is => 'Boolean',
            doc => 'report loci with zero methylation ratios',
            default => 1,
        },
        no_header => {
            is => 'Boolean',
            doc => 'do not put a header on the file',
            default => 0,
        },

# other options not exposed:
#   -u, --unique          process only unique mappings/pairs.
#   -p, --pair            process only properly paired mappings.
#   -q, --quiet           don't print progress on stderr.
#   -r, --remove-duplicate
#                         remove duplicated reads.
#   -t N, --trim-fillin=N
#                         trim N end-repairing fill-in nucleotides. [default: 2]
#   -g, --combine-CpG     combine CpG methylaion ratios on both strands.
#   -m FOLD, --min-depth=FOLD
#                         report loci with sequencing depth>=FOLD. [default: 1]

    ],
};

sub available_methratio_versions {
    return keys %METHRATIO_VERSIONS;
}

sub _reference_fasta {
    my ($self) = @_;

    my %mapping = (
        36 => 'NCBI-human-build36',
        37 => 'GRCh37-lite-build37',
    );

    my $reference = $self->reference;

    if (exists $mapping{$reference}) {
        my $reference_build = Genome::Model::Build::ReferenceSequence->get(
            name => $mapping{$reference}
        );
        $reference = $reference_build->cached_full_consensus_path('fa');
    }

    return $reference;
}

sub _generate_command_line {
    my ($self) = @_;

    my @cmd = (
        'python', $METHRATIO_VERSIONS{$self->version},
        '-o', File::Spec->join($self->output_directory, $self->output_file),
        '-d', $self->_reference_fasta
    );

    if ($self->output_zeros) {
        push @cmd, '-z';
    }

    if ($self->chromosome) {
        push @cmd, '-c', Genome::Sys->quote_for_shell($self->chromosome);
    }

    if ($self->no_header) {
        push @cmd, '-h';
    }

    push @cmd, $self->bam_file;

    return join ' ', @cmd;
}

sub execute {
    my ($self) = @_;

    if ($self->chromosome) {
        my $sanitized_chromosome = sanitize_string_for_filesystem(
            $self->chromosome);
        my $chromosome_output_dir = File::Spec->join(
            $self->output_directory, $sanitized_chromosome);
        Genome::Sys->create_directory($chromosome_output_dir);
        $self->output_directory($chromosome_output_dir);
    }

    unless (defined $self->version) {
        $self->error_message('methratio version %s not found.', $self->version);
        return 0;
    }

    my $return = Genome::Sys->shellcmd(
        cmd => $self->_generate_command_line,
        input_files => [$self->bam_file, $self->_reference_fasta],
    );

    unless ($return) {
        die $self->error_message('Failed to execute: returned %s', $return);
    }

    return 1;
}

1;
