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
            shell_args_position => 1,
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

    my $new_test_name = $self->new_test_name;
    for my $software_result ($self->software_results) {
        my $sr_id = $software_result->id;
        $self->status_message(
                "Setting test_name to \"$new_test_name\" on id:$sr_id");

        $software_result->test_name($new_test_name);
        if ($software_result->test_name eq $new_test_name) {
            $self->status_message("    Success! New lookup_hash is " .
                    $software_result->lookup_hash);
        } else {
            $self->status_message("    Failed! The test name is still \"" .
                    $software_result->test_name . "\"");
        }
    }

    return 1;
}

1;
