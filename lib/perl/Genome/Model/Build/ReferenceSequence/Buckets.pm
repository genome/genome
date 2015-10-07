package Genome::Model::Build::ReferenceSequence::Buckets;

use strict;
use warnings;

use Genome;
use Algorithm::Bucketizer;
use List::Util qw(max);

class Genome::Model::Build::ReferenceSequence::Buckets {
    is => 'Genome::SoftwareResult::StageableSimple',
    has_input => [
        reference_sequence_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'The reference to bucketize',
        },
    ],
    has_metric => [
        count => {
            is => 'Number',
            doc => 'The number of buckets for this reference',
        },
    ],
};

sub bucket_list {
    my $self = shift;

    my @buckets;
    for my $bucket (1..$self->count) {
        push @buckets, $self->bucket($bucket);
    }

    return \@buckets;
}

sub bucket {
    my $self = shift;
    my $bucket = shift;

    my $output_dir = $self->output_dir;

    my $bucket_file = File::Spec->join($output_dir, $bucket);
    my @items = Genome::Sys->read_file($bucket_file);
    chomp @items;

    return \@items;
}

sub _run {
    my $self = shift;

    my $chr_lengths = $self->_chromosomes_with_lengths;
    my $max_length = max( map $_->[1], @$chr_lengths );

    my $bucketizer = Algorithm::Bucketizer->new(bucketsize => $max_length);
    for my $chr_with_length (@$chr_lengths) {
        $bucketizer->add_item( @$chr_with_length );
    }

    #skip optimization to preserve order
    $bucketizer->optimize(maxrounds => 100_000);

    my @buckets = $bucketizer->buckets();
    $self->count(scalar(@buckets));

    for my $bucket (@buckets) {
        $self->_write_bucket($bucket);
    }

    return 1;
}

sub _chromosomes_with_lengths {
    my $self = shift;

    my $reference = $self->reference_sequence_build;
    my $seqdict = $reference->sequence_dictionary_path('sam');

    unless(-s $seqdict) {
        die $self->error_message('No sequence dictionary found for reference %s to create buckets', $reference->__display_name__);
    }

    my $parser = Genome::Utility::IO::SeparatedValueReader->create(
        separator => "\t",
        input => $seqdict,
        headers => ['tag', 'name', 'length'],
        allow_extra_columns => 1,
        ignore_lines_starting_with => '(?!@SQ)',
    );

    my @chr_lengths;

    while (my $line = $parser->next) {
        unless ($line->{tag} eq '@SQ') {
            Carp::confess 'parser error';
        }

        my ($sn, $chr) = split(':', $line->{name});
        unless ($sn eq 'SN') {
            die $self->error_message('Expected SN first in @SQ line but got %s', $sn);
        }

        my ($ln, $length) = split(':', $line->{length});
        unless($ln eq 'LN') {
            die $self->error_message('Expected LN second in @SQ line but got %s', $ln);
        }

        push @chr_lengths, [$chr, $length];
    }

    return \@chr_lengths;
}

sub _write_bucket {
    my $self = shift;
    my $bucket = shift;

    my $staging_dir = $self->temp_staging_directory;
    my $bucket_file = File::Spec->join($staging_dir, $bucket->serial);

    Genome::Sys->write_file($bucket_file, map "$_\n", $bucket->items);
}

1;
