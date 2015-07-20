package Genome::Model::Command::Admin::RemoveDiskAllocationsFromTestdb;

use strict;
use warnings;

use Genome;
use DBI;
use TestDbServer::CmdLine qw(search_databases get_template_by_id get_database_by_id);
use Genome::Config;

class Genome::Model::Command::Admin::RemoveDiskAllocationsFromTestdb {
    is => 'Command::V2',
    doc => 'Remove disk allocations that were created while running under a test database',
    has => [
        database_name => {
            is => 'Text',
            default_value => _get_default_database_name(),
            doc => 'test database name, derived from the ds_gmschema_server config value',
        },
        database_server => {
            is => 'Text',
            default_value => _get_default_database_server(),
            doc => 'database server name, derived from the ds_gmschema_server config value',
        
        },
        database_port => {
            is => 'Text',
            default_value => _get_default_database_port(),
            doc => 'database listening port, derived from the ds_gmschema_server config value',
        },
        template_name => {
            is => 'Text',
            default_value => _get_default_template_name(),
            doc => 'template database name, resolved from the database host in the ds_gmschema_server config value',
        },
    ],
};

sub _get_default_database_name {
    my $conn = _parse_database_connection_info_from_env_var();
    return $conn->{dbname};
}

sub _get_default_database_server {
    my $conn = _parse_database_connection_info_from_env_var();
    return $conn->{host};
}

sub _get_default_database_port {
    my $conn = _parse_database_connection_info_from_env_var();
    return( $conn->{port} // 5432 );
}

sub _get_default_template_name {
    my $test_db_name = __PACKAGE__->_get_default_database_name();
    return __PACKAGE__->get_template_name_from_database_name($test_db_name);
}

sub _parse_database_connection_info_from_env_var {
    my %connection;
    foreach my $key ( qw( dbname host port ) ) {
        no warnings 'uninitialized';
        ($connection{$key}) = $ENV{XGENOME_DS_GMSCHEMA_SERVER} =~ m/$key=(.*?)(?:;|$)/;
    }
    return \%connection;
}

sub execute {
    my $self = shift;

    unless ($self->is_running_in_test_env) {
        die $self->error_message('Must be run within a test environment created with genome-test-env');
    }

    my @allocations = $self->collect_newly_created_allocations();
    $self->report_allocations_to_delete(@allocations);
    $self->delete_allocations(@allocations);
    return 1;
}

sub is_running_in_test_env {
    my $self = shift;

}

my $sql_for_allocations = q(SELECT id FROM disk.allocation ORDER BY id);
sub collect_newly_created_allocations {
    my $self = shift;

    my $tmpl_dbh = $self->_dbh_for_template();
    my $db_dbh = $self->_dbh_for_database();

    my $tmpl_sth = $tmpl_dbh->prepare($sql_for_allocations);
    $tmpl_sth->execute();
    my $db_sth = $db_dbh->prepare($sql_for_allocations);
    $db_sth->execute();

    my @new_allocations_in_database;
    my($next_tmpl_allocation_id, $next_db_allocation_id)
    while(1) {
        unless (defined $next_tmpl_allocation_id) {
            ($next_tmpl_allocation_id) = $tmpl_sth->fetchrow_array();
        }
        unless (defined $next_db_allocation_id) {
            ($next_db_allocation_id) = $db_sth->fetchrow_array();
        }

        last unless defined($next_db_allocation_id);

        if ($next_tmpl_allocation_id eq $next_db_allocation_id) {
            # this allocation exists in both the template and test database, ignore it
            undef($next_tmpl_allocation_id);
            undef($next_db_allocation_id);

        } elsif ($next_tmpl_allocation_id lt $next_db_allocation_id) {
            # This allocation was deleted in the test database?!
            undef($next_tmpl_allocation_id);

        } else {
            # This allocation was created in the test database
            push @new_allocations_in_database, $next_db_allocation_id;
            undef($next_db_allocation_id);
        }
    }
    return @new_allocations_in_database;
}

sub _dbh_for_template {
    my $self = shift;
    return $self->_dbh_for('template');
}

sub _dbh_for_database {
    my $self = shift;
    return $self->_dbh_for('database');
}

sub _dbh_for {
    my($self, $db_or_tmpl) = @_;

    my $db_name = $db_or_tmpl eq 'database'
                    ? $self->database_name
                    : $self->template_name;
    my $db_host = $self->database_server;
    my $db_port = $self->database_port;

    return DBI->connect("dbi:Pg:dbname=${db_name};host=${db_host};port=${db_port}",
                        Genome::Config::get('ds_gmschema_login'),
                        Genome::Config::get('ds_gmschema_auth'),
                        { AutoCommit => 0, RaiseError => 1, PrintError => 0 });
}


sub report_allocations_to_delete {
    my $self = shift;

}

sub delete_allocations {
    my $self = shift;

}

sub get_template_name_for_database_name {
    my($self, $db_name) = @_;

    my @database_ids = search_databases(name => $db_name);
    unless (@database_ids == 1) {
        die $self->error_message('Expected 1 database named %s but got %d', $db_name, scalar(@database_ids));
    }

    my $database = get_database_by_id($databases[0]);
    my $tmpl = get_template_by_id($database->{template_id});
    return $tmpl->{name};
}

1;
