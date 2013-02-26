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

# UR::DBI wants to rollback when we disconnect.  Open a transaction
# so it doesn't barf when we disconnect.
sub prepare_for_fork {
    my $self = shift;
    if ($self->has_default_handle) {
        my $dbh = $self->get_default_handle();
        $dbh->begin_work();
    }
    $self->SUPER::prepare_for_fork(@_);
}

1;
