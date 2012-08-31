package Genome::Model::Tools::Annotate::Chromosome;

use strict;
use warnings;

use Genome;

my $DEFAULT_OUTPUT_FORMAT = 'gtf';
my $DEFAULT_VERSION = '54_36p_v2';
my $DEFAULT_ANNO_DB = 'NCBI-human.combined-annotation';

class Genome::Model::Tools::Annotate::Chromosome {
    is => ['Command'],
    has_input => [
        anno_db => {
            doc => 'The name of the annotation database to use. default_value='. $DEFAULT_ANNO_DB,
            default_value => $DEFAULT_ANNO_DB,
        },
        version => {
            doc => 'The version of the annotation database. default_value='. $DEFAULT_VERSION,
            default_value => $DEFAULT_VERSION,
        },
        chromosome => {
            #Should only be optional if run in parallel
            is_optional => 1,
            doc => 'The chromosome name to generate annotation for.',
        },
        output_format => {
            doc => 'The file format to output annotation in.',
            valid_values => ['gtf','gff','bed','gff3'],
            default_value => $DEFAULT_OUTPUT_FORMAT,
        },
        output_directory => {
            doc => 'The output directory where the annotation file will be dumped.',
        },
    ],
    has_output => [
        anno_file => {
            doc => 'The output file where annotation is written.  Do not define from command line',
            is_optional => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    Genome::Sys->create_directory($self->output_directory);
    $self->anno_file($self->output_directory .'/'. $self->anno_db .'_'. $self->version .'_'. $self->chromosome .'.'. $self->output_format);
    my $fh = Genome::Sys->open_file_for_writing($self->anno_file);
    unless ($fh) {
        die('Failed to open '. $self->output_format .' file: '. $self->anno_file);
    }
    my $ti = Genome::Model->get(name => $self->anno_db)->build_by_version($self->version)->transcript_iterator(chrom_name => $self->chromosome);
    my $format_string = $self->output_format .'_string';
    my %gene_strings;
    while (my $t = $ti->next) {
        my $gene = $t->gene;
        unless (defined($gene_strings{$gene->gene_id})) {
            my $gene_string = $gene->$format_string;
            if ($gene_string) {
                $gene_strings{$gene->gene_id} = $gene_string;
                print $fh $gene_string;
            } else {
                $gene_strings{$gene->gene_id} = 0;
            }
        }
        my $transcript_string = $t->$format_string;
        if ($transcript_string) {
            print $fh $transcript_string;
        }
        my @sub_structure = grep {$_->structure_type ne 'flank'} $t->ordered_sub_structures;
        for my $ss (@sub_structure){
            my $ss_string = $ss->$format_string;
            if ($ss_string) {
                print $fh $ss_string;
            }
        }
    }
    $fh->close;
    return 1;
}

1;
