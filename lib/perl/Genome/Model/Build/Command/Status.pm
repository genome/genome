package Genome::Model::Build::Command::Status;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::Status {
    is => 'Genome::Command::Base',
    doc => "prints status of non-succeeded builds and tallies all build statuses",
    has => [
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
            require_user_verify => 0,
            doc => 'Build(s) to check status. Resolved from command line via text string.',
            shell_args_position => 1,
        },
        hide_statuses => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            doc => 'Hide build details for statuses listed.',
        },
        show_event_info => {
            is => 'Boolean',
            default => 0,
            doc => 'BETA MAY NOT STAY - Show parent and oldest event info.',
        },
    ],
};

sub execute {
    my $self = shift;

    my @hide_statuses = $self->hide_statuses;

    my %status;
    my @builds = sort {$a->model_name cmp $b->model_name} $self->builds;
    my $model_name;
    for my $build (@builds) {
        my $build_status = $build->status;
        $status{$build_status}++;
        if (not grep { lc $_ eq lc $build_status } @hide_statuses) {
            if (!$model_name || $model_name ne $build->model_name) {
                $model_name = $build->model_name;
                $self->print_message("\nModel: ".$model_name);
                my $header = "BUILD_ID\tSTATUS";
                $header .= "\t\tPARENT_EVENT_ID\tPARENT_START_TIME\tOLDEST_CHILD_ID\tOLDEST_CHILD_START" if ($self->show_event_info);
                $self->print_message($header);
            }

            my $info = $build->id."\t$build_status";

            if ($self->show_event_info) {
                my $build_id = $build->id;
                my @running_events = Genome::Model::Event->get(build_id => $build_id, event_status => 'Running');
                @running_events = sort { $a->date_scheduled cmp $b->date_scheduled } @running_events;
                my ($parent_event) = grep { not defined $_->parent_event_id } @running_events;
                @running_events = grep { $_->id != $parent_event->id } @running_events;
                if ($parent_event) {
                    $info .= "\t\t" . $parent_event->id . "\t" . $parent_event->date_scheduled;
                }
                else {
                    $info .= "\t\t\t";
                }
                if (@running_events) {
                    $info .= "\t" . $running_events[0]->id . "\t" . $running_events[0]->date_scheduled;
                }
            }

            $self->print_message($info);
        }
    }

    my $total;
    for my $key (sort keys %status) {
        $total += $status{$key};
    }

    print "\n";
    for my $key (sort keys %status) {
        print "$key: $status{$key}\t";
    }
    print "Total: $total\n";

    return 1;
}

sub print_message {
    my $self = shift;
    my $msg = shift;
    print "$msg\n";
}

1;
