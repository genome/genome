package Genome::SoftwareResult::Command::SetTestName;

use strict;
use warnings;

use Genome;


class Genome::SoftwareResult::Command::SetTestName {
    is => 'Command::V2',

    has_input => [
        software_results => {
            is => 'Genome::SoftwareResult',
            is_many => 1,
            doc => 'software result to set name on'
        },
        new_test_name => {
            is => 'Text',
            doc => 'new test name to set',
        },
    ],
};

sub execute {
    my $self = shift;

    for my $software_result (@{$self->software_results}) {
        $software_result->set_test_name($self->new_test_name);
    }

    return 1;
}

1;
