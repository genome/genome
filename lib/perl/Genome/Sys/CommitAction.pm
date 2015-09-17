package Genome::Sys::CommitAction;

use Genome;

my $sequence = 0;
class Genome::Sys::CommitAction {
    has_optional => [
        on_sync     => { is => 'CODE', doc => 'Callback to run when data is saved to the data sources' },
        on_commit   => { is => 'CODE', doc => 'Callback to run when data sources are committed' },
        on_sync_fail => { is => 'CODE', doc => 'Callback to run when other on_sync callbacks have failed' },
        on_rollback => { is => 'CODE', doc => 'Callback to run when data sources are rolled-back' },
        data => { doc => 'A piece of data passed to the callbacks' },
    ],
    data_source => 'UR::DataSource::Default',
    id_generator => sub { ++$sequence },
};

my $rollback_mode = 'on_rollback';
foreach my $callback_data ( [ precommit => 'on_sync_fail' ], [ prerollback => 'on_rollback' ] ) {
    my($aspect, $mode) = @$callback_data;
    UR::Observer->register_callback(
        subject_class_name => 'UR::Context',
        subject_id => UR::Context->process->id,
        aspect => $aspect,
        callback => sub {
            $rollback_mode = $mode;
        }
    );
}

sub __save__ {
    my $self = shift;
    $self->_run('on_sync');
}

sub __commit__ {
    my $self = shift;
    $self->_run('on_commit');
    $self->delete();
}

sub __rollback__ {
    my $self = shift;
    $self->_run($rollback_mode);
    $self->delete() if $rollback_mode eq 'on_rollback';
}

sub _run {
    my($self, $which) = @_;
    if (my $cb = $self->$which) {
        $cb->($self->data);
    }
}

# These are never loaded, only created
sub __load__ {
    my($class, $bx, $headers) = @_;
    return($headers, []);
}

1;

__END__

=pod

=head1 NAME

Genome::Sys::CommitAction - Schedule code to run at commit time

=head1 SYNOPSIS

  Genome::Sys::CommitAction->create(
      on_sync => sub { write_temp_file($temp_file) },
      on_commit => sub { rename $temp_file, $permanent_name },
      on_sync_fail => sub { log_message("Commit failed, left temp file $temp_file") },
      on_rollback => sub { unlink($temp_file) },
  );

  # Later on...
  UR::Context->commit(); # on_sync then on_commit subs run here, when data sources are committed

=head1 DESCRIPTION

This class provides a mechanism to schedule some code to run later when
changed data is being saved to the database.  The callbacks are run only when
the base Process context is being committed or rolled-back, not when a
software transaction is committed/rolled-back.  Callbacks are only run once,
even if the Process context is committed or rolled-back multiple times.

There are two phases to committing saved data.  First, all changed data is
saved to the appropriate data source (colloquially called sync_database).
After all data is saved, each data source is asked to commit the changes.
If any data fails to save during sync_database, all data sources are asked
to rollback.

For these CommitAction objects, their C<on_sync> callback is run during
sync_database, C<on_commit> is run during commit.  C<on_rollback> is run when
the Process Context is rolled back.

The C<on_sync_fail> callback is only run when a subsequent C<on_sync> callback
throws an exception.  CommitActions are then considered live again, and
all their C<con_sync> callbacks may run again, if the original Context commit
exception is caught and it tries to commit again.  Note that C<on_sync_fail>
callbacks are _not_ run when the context is rolled back with
UR::Context->rollback.

The callbacks are passed a single argument, whatever is stored in the
CommitAction's 'data' attribute.

Also note that if the C<con_sync> or C<on_commit> callbacks change the program
state by creating, deleting or changing objects, those changes will not be
saved during the commit() the callbacks are running in; they will be saved
during the next commit.  The reason being that all the changes are collected
before any saving and committing is done.  Changes that happen during these
callbacks will have happened too late to make it into this list.  In
particular, this behavior differs from a C<precommit> callback on the Context
in that precommit callbacks run before the list of changes is collected.

Finally, these actions are implemented by creating ordinary UR Objects.  That
means if a CommitAction is created within a software transaction, it will
be deleted without running its callbacks if that transaction is rolled back.

=head1 SEE ALSO

L<UR::Context>

=cut
