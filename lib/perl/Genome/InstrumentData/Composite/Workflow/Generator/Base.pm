package Genome::InstrumentData::Composite::Workflow::Generator::Base;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Composite::Workflow::Generator::Base {
    is => 'UR::Singleton',
    is_abstract => 1,
};

my $counter = 0;
sub next_counter_value {
    return ++$counter;
}

sub _general_workflow_input_properties {
    my $class = shift;

    return qw(
        picard_version
        samtools_version
        bedtools_version
        force_fragment
        trimmer_name
        trimmer_version
        trimmer_params
        result_users
    );
}

1;
