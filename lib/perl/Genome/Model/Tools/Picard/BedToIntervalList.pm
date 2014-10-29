package Genome::Model::Tools::Picard::BedToIntervalList;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::BedToIntervalList {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input => {
            is  => 'String',
            doc => 'The input BED file',
        },
        output => {
            is  => 'String',
            doc => 'The output Picard Interval List',
        },
        sequence_dictionary => {
            is => 'String',
            doc => 'The sequence dictionary',
        },
    ],
};

sub help_brief {
    'Converts a BED file to an Picard Interval List.';
}

sub help_detail {
    return <<EOS

Documentation: http://broadinstitute.github.io/picard/command-line-overview.html#BedToIntervalList

Converts a BED file to an Picard Interval List.

EOS
}

sub execute {
    my $self = shift;
    
    my $jar_path = $self->picard_path .'/BedToIntervalList.jar';
    unless (-e $jar_path) {
        $self->error_message('Failed to find jar path: '. $jar_path);
        die($self->error_message);
    }
    my $cmd = $jar_path .' picard.util.BedToIntervalList INPUT='. $self->input .' OUTPUT='. $self->output .' SEQUENCE_DICTIONARY='. $self->sequence_dictionary;
    $self->run_java_vm(
        cmd => $cmd,
        input_files => [$self->input,$self->sequence_dictionary],
        skip_if_output_is_present => 0,
    );
    
    return 1;
}


1;
