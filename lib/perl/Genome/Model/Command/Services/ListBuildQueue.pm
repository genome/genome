package Genome::Model::Command::Services::ListBuildQueue;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Services::ListBuildQueue {
    is => 'Command',
    doc => 'List queued models for processing by cron.',
    has => [
        delimiter => {
            doc => 'delimiter used to separate model IDs per line',
            is => 'Text',
            default => "\n",
            is_optional => 1,
        },
        max => {
            doc => 'maximum number of scheduled builds per user',
            is => 'Integer',
            default => 50,
            is_optional => 1,
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
        {build_requested => '1', type_name => 'genotype microarray', -order_by => 'subject_id'},
        {build_requested => '1', -order_by => 'subject_id'},
    );

    my %count;
    while (my $iterator_params = shift @iterator_params) {
        my $models = Genome::Model->create_iterator(%{$iterator_params});
        while (my $model = $models->next) {
            my $username = $model->should_run_as() ? $model->run_as : Genome::Sys->username;

            unless (exists $count{$username}) {
                $count{$username} = scheduled_builds_for($username);
            }

            if ($count{$username} >= $self->max) {
                next;
            }

            $count{$username}++;
            $self->print($model->id);
        }
    }

    return 1;
}

sub scheduled_builds_for {
    my $username = shift;

    my $type = Genome::Model::Build->__meta__;

    my $table_name    = $type->table_name;
    my $run_by_column = $type->properties(property_name => 'run_by')->column_name;
    my $status_column = $type->properties(property_name => 'status')->column_name;

    my $sql_t = q/SELECT count(*) FROM %s WHERE %s = ? AND %s = 'Scheduled'/;
    my $sql = sprintf($sql_t, $table_name, $run_by_column, $status_column);

    my $dbh = $type->data_source->get_default_handle;
    my $sth = $dbh->prepare($sql);
    $sth->execute($username);
    my $row = $sth->fetchrow_arrayref();
    return $row->[0];
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
