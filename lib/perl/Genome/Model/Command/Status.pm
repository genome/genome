package Genome::Model::Command::Status;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Status {
    is => 'Genome::Command::Base',
    doc => "prints status of non-succeeded latest-builds and tallies latest-build statuses",
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            require_user_verify => 0,
            doc => 'Model(s) to check latest-build status. Resolved from command line via text string.',
            shell_args_position => 1,
        },
    ],
    has_optional => [
        summary_only => {
            is => 'Boolean',
            doc => "Only print the summary of the models' status.",
        },
    ],
};

sub execute {
    my $self = shift;

    my (%status, $total);
    for my $model ($self->models) {
        my $model_name = $model->name;
        my ($status, $build) = $model->status_with_build;
        my $build_id = ($build ? $build->id : 'N/A');
        $status{$status}++;
        $total++;
        print join("\t", $model_name, $build_id, $status) . "\n" if not $self->summary_only;
    }

    for my $key (sort keys %status) {
        print "$key: $status{$key}\t";
    }
    print "Total: $total\n";

    return 1;
}

1;

