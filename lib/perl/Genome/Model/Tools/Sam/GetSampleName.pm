package Genome::Model::Tools::Sam::GetSampleName;

use strict;
use warnings;

use above 'Genome';
use IPC::System::Simple qw(capture);
use Set::Scalar;

class Genome::Model::Tools::Sam::GetSampleName {
    is => 'Genome::Model::Tools::Sam',
    has_input => [
        bam_file => {
            is => 'Path',
            doc => 'The BAM file',
            shell_args_position => 1,
        },
    ],
    has_optional_output => [
        sample_name => {
            is => 'Text',
        },
    ],
};

sub execute {
    my $self = shift;

    die "bam_file could not be found" unless -e $self->bam_file;
    my $cmd = sprintf(q(%s view -H %s | grep '@RG'), $self->samtools_path, $self->bam_file);

    my @lines = capture($cmd);
    my $sample_names_set = Set::Scalar->new();
    for my $line (@lines) {
        $sample_names_set->insert(get_sample_name($line));
    }

    my @sample_names = $sample_names_set->members();
    my $num_sample_names = scalar(@sample_names);
    if ($num_sample_names == 0) {
        die sprintf("Found no sample names in %s", $self->bam_file);
    } elsif ($num_sample_names > 1) {
        die sprintf("Found multiple sample names in %s: %s", $self->bam_file, join(', ', @sample_names));
    } else {
        $self->sample_name($sample_names[0]);
        $self->status_message($self->sample_name);
    }
    return 1;
}

sub get_sample_name {
    my ($line) = @_;

    my %kv_pairs = get_key_value_pairs($line);
    return $kv_pairs{'SM'};
}

sub get_key_value_pairs {
    my $line = shift;

    my @pairs = split(/\t/, $line);
    my %result;
    for my $pair (@pairs) {
        my ($key, $value) = split(/:/, $pair);
        $result{$key} = $value;
    }
    return %result;
}

