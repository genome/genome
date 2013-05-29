package Genome::DataSource::PGTest;

use strict;
use warnings;
use Genome;
use Carp;

class Genome::DataSource::PGTest {
    is => 'Genome::DataSource::Main',
    has_constant => [
        server => { default_value => 'dbname=genome;host=gms-postgres' },
        login => { default_value => 'genome' },
        auth => { default_value => 'TGIlab' },
        owner => { default_value => 'public' },
    ],
};

sub _my_data_source_id {
    'Genome::DataSource::GMSchema';
}

sub _ds_tag {
    'Genome::DataSource::PGTest';
}

sub _dbi_connect_args {
    my $self = shift;
    
    my @connection = $self->SUPER::_dbi_connect_args(@_);
    
    $connection[3] = { AutoCommit => 1, RaiseError => 0 };

    return @connection;
}

# do not allow queries or connects when paused.

sub create_iterator_closure_for_rule {
  my $self = shift;
  
  $self->pause_db_queries_if_necessary();
  
  $self->SUPER::create_iterator_closure_for_rule(@_);
}

sub create_dbh {
  my $self = shift;
  $self->pause_db_queries_if_necessary();

  
  $self->SUPER::create_dbh(@_);
}

sub pause_db_queries_if_necessary {
    my $self = shift;
    return 1 unless $ENV{GENOME_DB_QUERY_PAUSE} and -e $ENV{GENOME_DB_QUERY_PAUSE};

    print "Database querying has been paused; disconnecting db handles.  Please wait until the query pause is released...\n";

    my @data_sources = UR::Context->all_objects_loaded('UR::DataSource::RDBMS');
    for my $ds (@data_sources) {
        $ds->disconnect_default_handle if $ds->has_default_handle;
    }

    while (1) {
        sleep 30;
        last unless $ENV{GENOME_DB_QUERY_PAUSE} and -e $ENV{GENOME_DB_QUERY_PAUSE};
    }

    print "Database updating has been resumed, continuing commit!\n";
    return 1;
}

sub _sync_database {
    my $self = shift;
    my $dbh = $self->get_default_handle();
   
    # do the actual sync in an eval, so we can set autocommit back
    my $sync_rv;
    eval {
        $dbh->begin_work(); 
        $sync_rv = $self->SUPER::_sync_database(@_);
    };
    
    # rethrow the exception if we got one
    if ($@) {
        die $@;
    }
    
    return $sync_rv;
}

# This is a no-op because commit is controlled by the filehandle trigger
# from our Oracle parent process.
sub commit {
    my $self = shift;
    return 1;
}

1;

