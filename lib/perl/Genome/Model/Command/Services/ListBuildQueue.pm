package Genome::Model::Command::Services::ListBuildQueue;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Services::ListBuildQueue {
    is => 'Command::V2',
    doc => 'List queued models for processing by cron.',
    has => [
        delimiter => {
            doc => 'delimiter used to separate model IDs per line',
            is => 'Text',
            default => "\n",
            is_optional => 1,
        },
        max => {
            doc => 'maximum number of scheduled builds per user (assumes run_as will be used)',
            is => 'Integer',
            default => 50,
            is_optional => 1,
        },
        running_max => {
            doc => 'maximum number of running builds per user (assumes run_as will be used)',
            default => 150,
        },
    ],
};

sub help_synopsis {
    return <<EOS;
genome model services list-build-queue | xargs -n 1 -P 5 /usr/local/bin/build-start
EOS
}

sub help_detail {
    return <<EOS;
List queued models for processing by cron.
EOS
}

# TODO is this needed?
# We had to debug a failure on the cron server where a downstream pipe crashed
# which meant that we got no mail due to having no SIGPIPE handler.
$SIG{PIPE} = sub {
    die "Cannot write to pipe (downstream pipe closed)";
};

sub execute {
    my $self = shift;

    my @iterator_params = (
        # prioritize genotype microarray over other builds because their
        # runtime is so short and are frequently prerequisite for other builds
        {build_requested => '1',  type_name     => 'genotype microarray', -order_by => 'subject_id'},
        {build_requested => '1', 'type_name !=' => 'genotype microarray', -order_by => 'subject_id'},
    );

    my @run_as = run_as_list();
    RUN_AS: while (my $run_as = shift @run_as) {
        my $scheduled_count = builds_for($run_as, 'Scheduled');

        # avoid extra DB queries below if we are just going to skip anyway
        if ($scheduled_count >= $self->max) {
            Genome::Logger->infof(
                "%s's scheduled count (%d) meets or exceeds specified maximum (%d)\n",
                $run_as,
                $scheduled_count,
                $self->max,
            );
            next RUN_AS;
        }

        my $max_running = $self->running_max;
        if(Genome::Sys->user_has_role($run_as, 'production')) {
            $max_running *= 5;
        }
        my $running_count = $scheduled_count + builds_for($run_as, 'Running');
        if($running_count > $max_running) {
            Genome::Logger->infof(
                "%s's running count (%d) meets or exceeds maximum (%d)\n",
                $run_as,
                $running_count,
                $max_running,
            );
            next RUN_AS;
        }

        ITER: for my $iterator_params (@iterator_params) {
            my $models = Genome::Model->create_iterator(%{$iterator_params}, run_as => $run_as);

            MODEL: while (my $model = $models->next) {
                if ($scheduled_count >= $self->max) {
                    Genome::Logger->infof(
                        "%s's scheduled count (%d) reached specified maximum (%d)\n",
                        $run_as,
                        $scheduled_count,
                        $self->max,
                    );
                    next RUN_AS;
                }
                if($running_count > $max_running) {
                    Genome::Logger->infof(
                        "%s's running count (%d) reached maximum (%d)\n",
                        $run_as,
                        $running_count,
                        $max_running,
                    );
                    next RUN_AS;
                }
                $scheduled_count++;
                $running_count++;
                $self->print($model->id);
            }
        }
    }

    return 1;
}

sub run_as_list {
    my $type = Genome::Model->__meta__;

    my $table_name    = $type->table_name;
    my $run_as_column = $type->properties(property_name => 'run_as')->column_name;
    my $build_requested_column = $type->properties(property_name => '_build_requested')->column_name;

    my $sql_t = q/SELECT distinct(%s) FROM %s WHERE %s IS TRUE/;
    my $sql = sprintf($sql_t, $run_as_column, $table_name, $build_requested_column);

    my $dbh = $type->data_source->get_default_handle;
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    my $rows = $sth->fetchall_arrayref();
    return map { $_->[0] } @{$rows};
}

sub builds_for {
    my $username = shift;
    my $status = shift;

    my $set = Genome::Model::Build->define_set(run_by => $username, status => $status);
    return $set->count();
}

sub sprint {
    my $self = shift;
    return join($self->delimiter, @_, '');
}

sub print {
    my $self = shift;
    print $self->sprint(@_);
}

1;
