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
        abandon_builds => {
            is => 'Boolean',
            default_value => 0,
        },
    ],
};

sub execute {
    my $self = shift;

    my $new_test_name = $self->new_test_name;

    my %builds;

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
        if ($self->abandon_builds) {
            my @builds = $software_result->builds;
            for my $build (@builds){
                if (!exists($builds{$build->id})) {
                    $builds{$build->id} = $build;
                }
            }
        }
    }
    if ($self->abandon_builds) {
        my @ids = keys %builds;
        my @builds = values %builds;
        $self->status_message('Abandoning builds: '. join(',',@ids));
        my $abandon_cmd = Genome::Model::Build::Command::Abandon->create(
            builds => \@builds,
            header_text => 'Build Abandoned - SR Test Name',
            body_text => 'Software results have new test name '. $self->new_test_name,
        );
        unless ($abandon_cmd) {
            $self->fatal_message('Unable to create abandon builds command!');
        }
        unless ($abandon_cmd->execute) {
            $self->fatal_messaage('Failed to execute abandon builds command!');
        }
    }

    return 1;
}

1;
