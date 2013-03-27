#eclark: Ignore table list should be reviewed.  I've seen classes for at few of them.

use strict;
use warnings;

package Genome::DataSource::GMSchemaOracle;

use Genome;

use Cwd;

class Genome::DataSource::GMSchemaOracle {
    is => ['UR::DataSource::Oracle'],
    type_name => 'genome datasource gmschema',
};

sub server {
    "dwrac";
}

sub login {
    "mguser";}

sub auth {
    "mguser_prd";}

sub owner {
    "MG";}

sub _my_data_source_id {
    "Genome::DataSource::GMSchema";
}


sub table_and_column_names_are_upper_case { 1; }

sub _sync_database {
    my $self = shift;

    $DB::single = 1;

    my $dbh = $self->get_default_handle;
    unless ($dbh->do("alter session set NLS_DATE_FORMAT = 'YYYY-MM-DD HH24:MI:SS'")
            and
            $dbh->do("alter session set NLS_TIMESTAMP_FORMAT = 'YYYY-MM-DD HH24:MI:SSXFF'"))
    {
        Carp::croak("Can't set date format: $DBI::errstr");
    }
    $self->SUPER::_sync_database(@_);
}

sub _resolve_class_name_for_table_name_fixups {
    my($self,@words) = @_;

    if ($words[0] eq 'Genome') {
        # Everything is already under Genome
        shift @words;
        if ($words[0] eq 'Model') {
            splice(@words, 1, 0, '::');  # Make it Model::Blah instead of ModelBlah
        }
    }
    return @words;
}

sub _lookup_class_for_table_name {
    my $self = shift;
    my $table_name = shift;

    # This currently depends on a bug in UR that will soon be fixed and will
    # then need to fall back onto Workflow::DataSource::InstanceSchema as a
    # consequence of our Postgres/Oracle syncing.
    my $class = $self->SUPER::_lookup_class_for_table_name($table_name);
    unless ($class) {
        my $gms = Genome::DataSource::GMSchema->get();
        $class = $gms->_lookup_class_for_table_name($table_name);
    }

    if (!$class && $ENV{GENOME_QUERY_POSTGRES}) {
        my %ora_map = Genome::DataSource::Main->oracle_to_postgres_table_mapping;
        my $mapped = $ora_map{lc($table_name)};
        print STDERR "OK falling back to $mapped\n";

        $class = $self->SUPER::_lookup_class_for_table_name($mapped);
    }
    return $class;
}


1;

