package Genome::DataSource::Meta;

# The datasource for metadata describing the tables, columns and foreign
# keys in the target datasource

use strict;
use warnings;

use version;

use UR;
use DBD::SQLite;

UR::Object::Type->define(
    class_name => 'Genome::DataSource::Meta',
    is => ['UR::DataSource::Meta'],
);

sub _dbi_connect_args {
    my $self = shift;

    my @connection = $self->SUPER::_dbi_connect_args(@_);

    if(version->parse($DBD::SQLite::VERSION) >= version->parse('1.38_01')) {
        my $connect_attr = $connection[3] ||= {};
        $connect_attr->{sqlite_use_immediate_transaction} = 0;
    }

    return @connection;
}

1;
