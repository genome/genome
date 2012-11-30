package Genome::SoftwareResult::Command::RemoveTestName;

use strict;
use warnings;

use Genome;


class Genome::SoftwareResult::Command::RemoveTestName {
    is => 'Command::V2',

    has_input => [
        software_results => {
            is => 'Genome::SoftwareResult',
            doc => 'software result to set name on',
            is_many => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    foreach my $sr ($self->software_results) {
        unless ($sr->remove_test_name) {
            $self->error_message("Failed to remove test_name for ".$sr->id);
        }
    }

    return 1;
}

1;
