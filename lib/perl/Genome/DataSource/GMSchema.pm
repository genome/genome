package Genome::DataSource::GMSchema;

use strict;
use warnings;
use Genome;

use File::Basename;
use File::Spec;
use Genome::DataSource::CommonRDBMS qw(log_error log_commit_time);

class Genome::DataSource::GMSchema {
    is => [$ENV{GENOME_DS_GMSCHEMA_TYPE}, 'Genome::DataSource::CommonRDBMS'],
    has_constant => [
        server  => { default_value  => $ENV{GENOME_DS_GMSCHEMA_SERVER} },
        login   => { default_value  => $ENV{GENOME_DS_GMSCHEMA_LOGIN} },
        auth    => { default_value  => $ENV{GENOME_DS_GMSCHEMA_AUTH} },
        owner   => { default_value => undef, is_optional => 1 },
    ],
};

if ($ENV{GENOME_TEST_FILLDB}) {
    # GENOME_TEST_FILLDB should be a complete DBI connect string, with user and password, like:
    # dbi:Pg:dbname=testdb;host=computername;port=5434;user=testuser;password=testpasswd
    Genome::DataSource::GMSchema->alternate_db_dsn( $ENV{GENOME_TEST_FILLDB} );
}

sub _ignore_table {
    my($self, $table_name) = @_;
    return 1 if $table_name =~ m/^sqitch/;
    return 1 if $table_name =~ m/^pg_catalog/;
    return 1 if $table_name =~ m/^information_schema/;
    return 1 if $table_name eq 'public.pg_stat_statements';
    return 1 if $table_name eq 'public.appmanager';
    return 0;
}


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
