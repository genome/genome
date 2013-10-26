package Genome::DataSource::GMSchema;

use strict;
use warnings;
use Genome;

use Genome::DataSource::CommonRDBMS qw(log_error log_commit_time);

class Genome::DataSource::GMSchema {
    is => ['UR::DataSource::Pg', 'Genome::DataSource::CommonRDBMS'],
    has_constant => [
        server => { default_value => 'dbname=genome' },
        login => { default_value => 'genome' },
        auth => { default_value => 'changeme' },
        owner => { default_value => 'public' },
    ],
};

sub table_and_column_names_are_upper_case { 0 }

sub _dbi_connect_args {
    my $self = shift;

    my @connection = $self->SUPER::_dbi_connect_args(@_);

    my $connect_attr = $connection[3] ||= {};
    $connect_attr->{AutoCommit} = 1;
    $connect_attr->{RaiseError} = 0;

    return @connection;
}

1;
