package Genome::Model::Tools::Picard::BedToIntervalList;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::BedToIntervalList {
    is  => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        input => {
            is  => 'String',
            doc => 'The input BED file',
            picard_param_name => 'INPUT',
        },
        output => {
            is  => 'String',
            doc => 'The output Picard Interval List',
            picard_param_name => 'OUTPUT',
        },
        sequence_dictionary => {
            is => 'String',
            doc => 'The sequence dictionary',
            picard_param_name => 'SEQUENCE_DICTIONARY',
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

sub _jar_name { 'BedToIntervalList.jar' }
sub _java_class_name { 'picard.util.BedToIntervalList' }

sub _shellcmd_extra_params {
    my $self = shift;
    return (
        input_files => [$self->input, $self->sequence_dictionary],
        skip_if_output_is_present => 0,
        );
}

sub _validate_params {
    my $self = shift;
    $self->enforce_minimum_version('1.120');
}

1;
