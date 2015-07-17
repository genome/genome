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

sub collect_newly_created_allocations {
    my $self = shift;

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
