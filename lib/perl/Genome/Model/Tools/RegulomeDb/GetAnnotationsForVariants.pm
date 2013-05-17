package Genome::Model::Tools::RegulomeDb::GetAnnotationsForVariants;

use strict;
use warnings;
use Genome;
use WWW::Mechanize;

class Genome::Model::Tools::RegulomeDb::GetAnnotationsForVariants {
    is => 'Command::V2',
    has => [
        variant_list => {
            is => 'String',
            doc => 'Path to variant list in bed format',
        },
        output_file => {
            is => 'String',
            doc => 'Path of output file',
        },
        format => {
            is => 'String',
            doc => 'format of regulomeDb annotation to get',
            default => 'full',
            valid_values => ['full', 'bed', 'gff'],
        },
    ],
};

sub execute {
    my $self = shift;

    my $content = Genome::Sys->read_file($self->variant_list);
    my $output = Genome::Model::Tools::RegulomeDb::fetch_large_annotation(
       $self->format, $content 
    );
    Genome::Sys->write_file($self->output_file, $output);
    return 1;
}

1;

