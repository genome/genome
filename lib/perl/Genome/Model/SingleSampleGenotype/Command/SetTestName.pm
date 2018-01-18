package Genome::Model::SingleSampleGenotype::Command::SetTestName;

use strict;
use warnings;

use Genome;

use File::Spec;

class Genome::Model::SingleSampleGenotype::Command::SetTestName {
    is        => 'Command::V2',
    doc       => 'Set a test name for all builds resolved',
    has_input => {
        builds => {
            is      => "Genome::Model::Build::SingleSampleGenotype",
            is_many => 1,
        },
        new_test_name => {
            is  => 'String',
            doc => 'The new test name string to set on all results'
        },
        abandon_builds => {
            is => 'Boolean',
            doc => 'attempt to abandon all of the builds used by these results',
            default_value => 0,
        },
    },
};

sub execute {
    my $self = shift;

    my @results_to_set_test_name;
    for my $build ( $self->builds ) {
        my @disk_usage_results = $build->disk_usage_results;
        push @results_to_set_test_name, @disk_usage_results;
    }
    my $set_test_name_cmd = Genome::SoftwareResult::Command::SetTestName->create(
        software_results => \@results_to_set_test_name,
        new_test_name => $self->new_test_name,
        abandon_builds => $self->abandon_builds,
    );
    unless ($set_test_name_cmd) {
        $self->fatal_message('Failed to create set test name command!');
    }
    unless ($set_test_name_cmd->execute) {
        $self->fatal_message('Unable to execute set test name command for all build software results!');
    }
    return 1;
}


1;
