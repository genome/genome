package Genome::Model::Tools::Mutect::ParallelWrapper;

use Data::Dumper;
use Genome;
use Carp qw/confess/;
use strict;
use warnings;

class Genome::Model::Tools::Mutect::ParallelWrapper {
    is => ['Genome::Model::Tools::Mutect'],
    doc => "Run the Mutect tool.",
    has_input => [
        fasta_object => {
            is => 'Genome::File::Fasta',
        },
        chunk_num => {
            is => 'Integer',
            doc => 'which chunk this is',
        },
        total_chunks => {
            is => 'Number',
            doc => 'total number of chunks to divide the genome into',
        },
        basename => {
            is => 'String',
            doc => 'basename for all files',
            is_optional => 1,
        },
    ],
    has_output => [
        output_file => {
            is => 'String',
            doc => 'Native output format file',
            is_optional => 1,
        },
    ],
};

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
This command exists in order to facilitate parallelization of
music's clinical correlation command. The output file will
be named \$output_directory/\$variant_matrix.results.glm.tsv.
This allows us to split a variant matrix into many parts and
process them in parallel easily with Workflow.
EOS
}

sub execute {
    my $self = shift;

    my $intervals = $self->convert_chunks_to_mutect_regions($self->fasta_object->divide_into_chunks($self->total_chunks, $self->chunk_num));
    $self->intervals($intervals);
    $self->output_file($self->basename . "_" . $self->chunk_num . ".out");
    $self->vcf($self->basename . "_" . $self->chunk_num . ".vcf");
    $self->coverage_file($self->basename . "_" . $self->chunk_num . ".wig");
    my $sub = $self->super_can('_execute_body');
    return $sub->($self, @_);
};

sub convert_chunks_to_mutect_regions {
    my ($self, @chunks) = @_;
    #this thing is a list of array refs 
        my @intervals;
        for my $interval (@chunks) {
            my ($chr, $start, $stop) = @$interval;
            push @intervals, "$chr:$start-$stop";
        }
    return \@intervals;
}

1;
