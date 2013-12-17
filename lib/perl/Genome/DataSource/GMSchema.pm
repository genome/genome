package Genome::DataSource::GMSchema;

use strict;
use warnings;
use Genome;

use Genome::DataSource::CommonRDBMS qw(log_error log_commit_time);

class Genome::DataSource::GMSchema {
    is => [$ENV{GENOME_DS_GMSCHEMA_TYPE}, 'Genome::DataSource::CommonRDBMS'],
    has_constant => [
        server  => { default_value  => $ENV{GENOME_DS_GMSCHEMA_SERVER} },
        login   => { default_value  => $ENV{GENOME_DS_GMSCHEMA_LOGIN} },
        auth    => { default_value  => $ENV{GENOME_DS_GMSCHEMA_AUTH} },
        owner   => { default_value  => $ENV{GENOME_DS_GMSCHEMA_OWNER} },
    ],
};

sub table_and_column_names_are_upper_case { 0 }

sub _dbi_connect_args {
    my $self = shift;

    my @connection = $self->SUPER::_dbi_connect_args(@_);

    my $connect_attr = $connection[3] ||= {};
    $connect_attr->{AutoCommit} = 0;  # this is set to true on master for unknown reasons
    $connect_attr->{RaiseError} = 0;

    return @connection;
}

1;
