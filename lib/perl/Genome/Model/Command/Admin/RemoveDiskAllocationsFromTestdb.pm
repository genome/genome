package Genome::Model::Command::Admin::RemoveDiskAllocationsFromTestdb;

use strict;
use warnings;

use Genome;
use DBI;

$DB::single=1;
class Genome::Model::Command::Admin::RemoveDiskAllocationsFromTestdb {
    is => 'Command::V2',
    doc => 'Remove disk allocations that were created while running under a test database',
    has => [
        database_name => {
            is => 'Text',
            default_value => _parse_database_name_from_env_var(),
            doc => 'test database name, derived from XGENOME_DS_GMSCHEMA_SERVER env var',
        },
        database_server => {
            is => 'Text',
            default_value => _parse_database_server_from_env_var(),
            doc => 'database server name, derived from XGENOME_DS_GMSCHEMA_SERVER env var',
        
        },
        database_port => {
            is => 'Text',
            default_value => _parse_database_port_from_env_var(),
            doc => 'database listening port, derived from XGENOME_DS_GMSCHEMA_SERVER env var',
        },
        template_name => {
            is => 'Text',
            default_value => _resolve_template_name_from_env_var(),
            doc => 'template database name, resolved from the database host in env XGENOME_DS_GMSCHEMA_SERVER env var',
        },
    ],
};

sub _parse_database_name_from_env_var {
    my $conn = _parse_database_connection_info_from_env_var();
    return $conn->{dbname};
}

sub _parse_database_server_from_env_var {
    my $conn = _parse_database_connection_info_from_env_var();
    return $conn->{host};
}

sub _parse_database_port_from_env_var {
    my $conn = _parse_database_connection_info_from_env_var();
    return $conn->{port};
}

sub _resolve_template_name_from_env_var {
    return 1;
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

1;
