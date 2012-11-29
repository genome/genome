package Genome::SoftwareResult::Command::SetTestName;

use strict;
use warnings;

use Genome;


class Genome::SoftwareResult::Command::SetTestName {
    is => 'Command::V2',

    has_input => [
        software_result => {
            is => 'Genome::Softwareresult',
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

    $self->software_result->set_test_name($self->new_test_name);

    return 1;
}

1;
