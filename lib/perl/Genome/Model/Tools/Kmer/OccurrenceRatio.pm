package Genome::Model::Tools::Kmer::OccurrenceRatio;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Kmer::OccurrenceRatio {
    is => 'Genome::Model::Tools::Kmer',
    has => [
        output_file => {
            is => 'Text',
            doc => 'The output file to print ratios',
        },
        output_type => {
            is => 'Text',
            doc => 'Specify what to output by giving at least one of the four keywords: unique, nonunique, nonuniquemulti, relative, and total.',
            default_value => 'nonunique unique relative',
        },
        index_name => {
            is => 'Text',
        },
        minimum_mer_size => {
            is => 'Integer',
            default_value => 1,
        },
        maximum_mer_size => {
            is => 'Integer',
            default_value => 32,
        },
        scan => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'When set the suffix arrary/tables are not loaded into memory but rather scanned reducing the memory requirements',
        },
    ],
};

sub execute {
    my $self = shift;

    my $gt_path = $self->genometools_path;
    my $options = '-output '. $self->output_type;
    if ($self->minimum_mer_size) {
        $options .= ' -minmersize '. $self->minimum_mer_size;
    }
    if ($self->maximum_mer_size) {
        $options .= ' -maxmersize '. $self->maximum_mer_size;
    }
    if ($self->scan) {
        $options .= ' -scan';
    }
    my $cmd = $gt_path .' tallymer occratio '. $options .' -esa '. $self->index_name .' > '. $self->output_file;
    my @output_files = ($self->output_file);
    Genome::Sys->shellcmd(
        cmd => $cmd,
        output_files => \@output_files,
    );
    return 1;
}
