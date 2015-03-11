package Genome::Sys::CommitAction;

use Genome;

my $sequence = 0;
class Genome::Sys::CommitAction {
    has_optional => [
        on_sync     => { is => 'CODE', doc => 'Callback to run when data is saved to the data sources' },
        on_commit   => { is => 'CODE', doc => 'Callback to run when data sources are committed' },
        on_sync_fail => { is => 'CODE', doc => 'Callback to run when other on_sync callbacks have failed' },
    ],
    data_source => 'UR::DataSource::Default',
    id_generator => sub { ++$sequence },
};

sub __save__ {
    my $self = shift;
    $self->_run('on_sync');
}

sub __commit__ {
    my $self = shift;
    $self->_run('on_commit');
}

sub __rollback__ {
    my $self = shift;
    $self->_run('on_sync_fail');
}

sub _run {
    my($self, $which) = @_;
    if (my $cb = $self->$which) {
        $cb->();
    }
}

1;

__END__

=pod

=head1 NAME

Genome::Sys::CommitAction - Schedule code to run at commit time

=head1 SYNOPSIS

  Genome::Sys::CommitAction->create(
      on_sync => sub { write_temp_file($temp_file) },
      on_commit => sub { rename $temp_file, $permanent_name }
      on_sync_fail => sub { log_message("Commit failed, left temp file $temp_file") }
  );

  # Later on...
  UR::Context->commit(); # on_commit sub runs here, when data sources are committed

=head1 DESCRIPTION

This class provides a mechanism to schedule some code to run later when
changed data is being saved to the database.  The callbacks are only when
the base, Process context is being committed, not when a software transaction
is committed.  Callbacks are only run once, even if the Process context is
committed multiple times.

There are two phases to committing saved data.  First, all changed data is
saved to the appropriate data source (colloquially called sync_database).
After all data is saved, each data soruce is asked to commit the changes.
If any data fails to save during sync_database, all data sources are asked
to rollback.

For these CommitAction objects, their C<on_sync> callback is run during
sync_database, C<on_commit> is run during commit.

The C<on_sync_fail> callback is only run when a subsequent C<on_sync> callback
throws an exception.  CommitActions are then considered live again, and
all their C<con_sync> callbacks may run again, if the original Context commit
exception is caught and it tries to commit again.  Note that C<on_sync_fail>
callbacks are _not_ run when the context is rolled back with
UR::Context->rollback

=head1 SEE ALSO

L<UR::Context>

=cut
