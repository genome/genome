package Genome::DataSource::CommonRDBMS;

use DBD::Pg v2.19.3;
use Genome;
use UR::DataSource::Pg;


class Genome::DataSource::CommonRDBMS {
    doc => 'Mixin class to implement pausing access to the database',
    valid_signals => [qw( precommit precreate_handle query sequence_nextval )]
};

my $query_pause = _make_db_pause_function('query_pause_sentry_file_path');
foreach my $signal ( qw( query precreate_handle sequence_nextval ) ){
    UR::Observer->register_callback(
        subject_class_name => __PACKAGE__,
        aspect => $signal,
        callback => $query_pause,
    );
}

my $commit_pause = _make_db_pause_function('commit_pause_sentry_file_path');
UR::Observer->register_callback(
    subject_class_name => 'UR::Context',
    subject_id => UR::Context->current->id,
    aspect => 'precommit',
    callback => sub {
        __PACKAGE__->$commit_pause;
    }
);

sub query_pause_sentry_file_path {
    Genome::Config::get('db_query_pause');
}

sub commit_pause_sentry_file_path {
    Genome::Config::get('db_pause');
}

sub pause_sleep_length_seconds { 30 }

sub _make_db_pause_function {
    my $pause_check_method = shift;

    return sub {
        my $self = shift;

        my $sentry = $self->$pause_check_method;
        return 1 unless $sentry and -e $sentry;

        print STDERR "Database querying has been paused; disconnecting db handles.  "
                    . "Please wait until the query pause is released...\n";

        my @data_sources = UR::Context->all_objects_loaded(__PACKAGE__);
        for my $ds (@data_sources) {
            $ds->disconnect_default_handle if $ds->has_default_handle;
        }

        while (1) {
            sleep pause_sleep_length_seconds;
            last unless $sentry and -e $sentry;
        }

        print STDERR "Database updating has been resumed, continuing query!\n";
        return 1;
    };
}


1;
