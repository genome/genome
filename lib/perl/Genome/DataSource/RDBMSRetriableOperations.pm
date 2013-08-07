package Genome::DataSource::RDBMSRetriableOperations;

use Genome;

# A mixin class that provides methods to retry queries and syncs
#
# Consumers should provide should_retry_operation_after_error().
# It's passed the SQL that generated the error and the DBI error string.
# It should return true if the operation generating that error should be
# retried.

class Genome::DataSource::RDBMSRetriableOperations { };


# The guts of the thing.  Consumers that want a base-datasource method to
# be retriable should override the method to call this instead, and pass
# a code ref to perform the retriable action

sub _retriable_operation {
    my $self = shift;
    my $code = shift;

    my $db_retry_time = 1; # seconds

    $self->_make_retriable_operation_observer();

    RETRY_LOOP:
    for( my $db_retry_time = 1; $db_retry_time < 3600; $db_retry_time *= 2 ) {
        my @rv = eval { $code->(); };

        if ($@) {
            if ($@ =~ m/DB_RETRY/) {
                $self->debug_message("Disconnecting and sleeping for $db_retry_time seconds...\n");
                $self->disconnect_default_handle;
                sleep $db_retry_time;
                next RETRY_LOOP;
            }
            Carp::croak($@);  # re-throw other exceptions
        }
        return $self->context_return(@rv);
    }
    die "Maximum database retries reached";
}


{
    my @retry_observers;
    sub _make_retriable_operation_observer {
        my $self = shift;
        unless (@retry_observers) {
            @retry_observers = map {
                $self->add_observer(
                    aspect => $_,
                    priority => 99999, # Super low priority to fire last
                    callback => \&_db_retry_observer,
                );
            }
            qw(query_failed commit_failed);
        }
    }
}

# Default is to not retry
sub should_retry_operation_after_error {
    my($self, $sql, $dbi_errstr) = @_;
    return 0;
}


# The callback for the retry observer
sub _db_retry_observer {
    my($self, $aspect, $db_operation, $sql, $dbi_errstr) = @_;

    $self->error_message("SQL failed during $db_operation\nerror: $dbi_errstr\nsql: $sql");

    die "DB_RETRY" if $self->should_retry_operation_after_error($sql, $dbi_errstr);

    # just fall off the end here...
    # Code triggering the observer will throw an exception
}

1;
