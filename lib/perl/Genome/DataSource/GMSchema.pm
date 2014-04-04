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
        owner   => { default_value  => $ENV{GENOME_DS_GMSCHEMA_OWNER} },
    ],
};

#BEGIN {
#    my($name, $path, $suffix) = File::Basename::fileparse($0, '.pm','.t','.pl');
#    $path .= '.testdb';
#    my $testdb_dirname = File::Spec->catdir($path, $name) . '/';
#    my $testdb_dirname = $0 . '.testdb/';
    
    if ($ENV{GENOME_TEST_FILLDB}) {
#        $ENV{UR_TEST_FILLDB} = "dbi:SQLite:dbname=$testdb_dirname";
#        $ENV{UR_TEST_FILLDB} = "dbi:Pg:dbname=genometest;host=10.0.108.127;user=genome;password=genome";
#        Genome::DataSource::GMSchema->alternate_db_dsn( "dbi:Pg:dbname=genometest;host=10.0.108.127;user=genome;password=genome");
        Genome::DataSource::GMSchema->alternate_db_dsn( "dbi:Pg:dbname=genometest;host=linus146;port=5434;user=genome;password=genome");
    }
    if ($ENV{GENOME_TEST_USEDB}) {
        $ENV{GENOME_DS_GMSCHEMA_TYPE} = 'UR::DataSource::Pg';
        $ENV{GENOME_DS_GMSCHEMA_SERVER} = 'dbname=genometest;host=10.0.108.127';
        $ENV{GENOME_DS_GMSCHEMA_LOGIN} = 'genome';
        $ENV{GENOME_DS_GMSCHEMA_AUTH} = 'genome';
#        $ENV{GENOME_DS_GMSCHEMA_TYPE} = 'UR::DataSource::SQLite';
#        $ENV{GENOME_DS_GMSCHEMA_SERVER} = $testdb_dirname;
    }
#}


{
    my %table_name_to_class_name_map = (
        'config.analysis_menu_item' => 'Genome::Config::AnalysisMenu::Item',
    );

    sub _lookup_class_for_table_name {
        my($self, $table_name) = @_;

        my $class_name = $table_name_to_class_name_map{$table_name};
        
        if ($class_name) {
            $class_name->class;
        } else {
            $class_name = $self->SUPER::_lookup_class_for_table_name($table_name);
        }
        return $class_name;
    }
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
