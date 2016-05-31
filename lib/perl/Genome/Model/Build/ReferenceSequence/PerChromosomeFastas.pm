package Genome::Model::Build::ReferenceSequence::PerChromosomeFastas;

use strict;
use warnings;

use Genome;
use File::Spec;
use Cwd;

class Genome::Model::Build::ReferenceSequence::PerChromosomeFastas {
    is => 'Genome::SoftwareResult::StageableSimple',
    has_input => [
        reference_sequence_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'The reference to generate per chromosome FASTA files for.',
        },
    ],
};

sub _run {
    my $self = shift;

    my $full_consensus = $self->reference_sequence_build->full_consensus_path('fa');
    my $cwd = getcwd();
    chdir($self->temp_staging_directory);
    my $cmd = 'awk \'BEGIN { CHROM="" } { if ($1~"^>") CHROM=substr($1,2); print $0 > CHROM".fa" }\' '. $full_consensus;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$full_consensus],
    );
    chdir($cwd);

    return 1;
}

sub fasta_directory {
    my $self = shift;

    my ($allocation) = $self->disk_allocations;
    return $allocation->absolute_path;
}

sub fasta_list {
    my $self = shift;

    my @fasta_files = glob(File::Spec->join($self->fasta_directory,'*.fa'));
    return \@fasta_files;
}


1;
