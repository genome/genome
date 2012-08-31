package Genome::Model::Tools::Abyss::Parallel;

use strict;
use warnings;

use Genome;
use File::Path qw/make_path/;
use File::pushd;

class Genome::Model::Tools::Abyss::Parallel {
    is => 'Genome::Model::Tools::Abyss',
    has_input => [
        kmer_size => {
            is => 'Number',
            doc => 'k-mer size for assembly',
        },
        min_pairs => {
            is => 'Number',
            doc => 'The minimum number of pairs needed to consider joining two contigs',
            default => 10,
        },
        min_sequence_identity => {
            is => 'Float',
            doc => 'Minimum sequence identity for PopBubbles and PathConsensus (real number in [0.0,1.0])',
            is_optional => 1
        },
        num_jobs => {
            is => 'Number',
            doc => 'The number of jobs to run in parallel',
            default => 8,
        },
        fastq_a => {
            is => 'Text',
            doc => 'fastq file a',
        },
        fastq_b => {
            is => 'Text',
            doc => 'fastq file b',
        },
        single_end_inputs => {
            is => 'Text',
            doc => 'Fasta or fastq files of single-end reads in addition to the paired inputs (space delimited, quote it!)',
            is_optional => 1,
        },
        name => {
            is => 'Text',
            doc => 'Name to prepend to output files',
            default => 'abyss',
        },
        job_queue => {
            is => 'Text',
            doc => 'The job queue to schedule the work in.',
            default => 'apipe',
        },
        min_coverage => {
            is => 'Text',
            doc => 'Remove contigs with mean k-mer coverage less than this value',
            is_optional => 1,
        },
        output_directory => {
            is => 'Text',
            doc => 'Directory to write output to',
        },
    ]
};

sub abyss_pe_binary {
    return shift->bindir . "/abyss-pe";
}

# TODO: abstract job submission, it's bad to hard code bsub all over
sub mpirun_cmd {
    my ($self, $output_directory) = @_;
    my $rusage = 'span[ptile=1]';
    my $log_file = "$output_directory/abyss_parallel.log";
    my $num_jobs = $self->num_jobs;
    my $job_queue = $self->job_queue;
    return "bsub -K -oo $log_file -n $num_jobs -a openmpi -q $job_queue -R '$rusage' mpirun.lsf -x PATH";
}

sub parse_kmer_range {
    my ($self, $range) = @_;

    my @kmer_sizes;
    if (my ($start, $end, $junk, $step) = $range =~ /^(\d+)\.\.(\d+)(:(-{0,1}\d+)|)/) {
        die "invalid kmer size range '$range', end=$end < start=$start. abort." if $end < $start;
        $step = 1 unless defined $step;
        die "invalid kmer size range '$range', step=$step <= 0, abort." if $step <= 0;
        for (my $i = $start; $i <= $end; $i += $step) {
            push(@kmer_sizes, $i);
        }
    } else {
        die "Invalid number '$range'" unless $range =~ /^\d+$/;
        push(@kmer_sizes, $range);
    }
    return @kmer_sizes;
}

sub get_kmer_sizes {
    my ($self, $kmer_size) = @_;
    my @kmer_sizes;
    my @values = split(',', $kmer_size);
    for my $v (@values) {
        push(@kmer_sizes, $self->parse_kmer_range($v));
    }
    return @kmer_sizes;
}

sub make_cmdline_args {
    my ($self, $kmer_size, $output_dir) = @_;

    my @cmd = (
        "v=-v",
        "k=" . $kmer_size,
        "n=" . $self->min_pairs,
        "np=" . $self->num_jobs,
        "in='" . $self->fastq_a . " " . $self->fastq_b . "'", 
        "name=".$self->name,
        'mpirun="' . $self->mpirun_cmd($output_dir).'"',
        );

    push(@cmd, "c=".$self->min_coverage) if $self->min_coverage;
    push(@cmd, "se='".$self->single_end_inputs."'") if $self->single_end_inputs;
    push(@cmd, "p=".$self->min_sequence_identity."") if $self->min_sequence_identity;

    return @cmd;
}

sub execute {
    my $self = shift;
    my $bindir = $self->bindir;

    die "Input fastq_a file '" . $self->fastq_a . " does not exist " unless -e $self->fastq_a;
    die "Input fastq_b file '" . $self->fastq_b . " does not exist " unless -e $self->fastq_b;

    for my $kmer_size ($self->get_kmer_sizes($self->kmer_size)) {
        my $output_dir = $self->output_directory . '/k' . $kmer_size;
        my $main_log_file = "$output_dir/abyss.log";
        my @cmd = (
            $self->abyss_pe_binary,
            $self->make_cmdline_args($kmer_size, $output_dir),
            " > $main_log_file 2>&1",
            );

        make_path($output_dir);
        pushd($output_dir);
        chdir($output_dir);
        local $ENV{PATH} = $self->bindir . ":$ENV{PATH}";
        return unless Genome::Sys->shellcmd( cmd => join(' ', @cmd),);
    }

    return 1;
}

1;
