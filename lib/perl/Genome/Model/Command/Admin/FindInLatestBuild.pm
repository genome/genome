package Genome::Model::Command::Admin::FindInLatestBuild;

class Genome::Model::Command::Admin::FindInLatestBuild {
    is => 'Command::V2',
    doc => 'Find models whose latest builds have <search text> in their .err or .out files.',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            shell_args_position => 1,
            require_user_verify => 0,
        },
        search_string => {
            is => 'Text',
            doc => 'Text you wish to search for',
        }
    ],
};

use strict;
use warnings;
use Genome;


sub execute {
    my ($self) = @_;

    my $search_string = $self->search_string;
    my $CMD = "grep -r -c -l -m 1 --include '*.err' --include '*.out' '$search_string'";

    my @models = $self->models;
    printf STDERR "Found %d models.\n", scalar(@models);

    my @affected_models;
    for my $model (@models) {
        my $latest_build = $model->latest_build;
        next unless $latest_build;

        my $data_directory = $latest_build->data_directory;
        next unless $data_directory;

        my @affected_logs = `$CMD $data_directory`;
        if (scalar(@affected_logs) > 0) {
            push @affected_models, $model;
            print STDOUT $model->id . "\n";
        } else {
            print STDERR "Unaffected " . $model->id . "\n";
        }
    }
    printf STDERR "%d/%d Models were affected.\n", scalar(@affected_models), scalar(@models);
}
