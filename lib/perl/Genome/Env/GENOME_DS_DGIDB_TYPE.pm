package Genome::Env::GENOME_DS_DGIDB_TYPE;

use Carp;

unless ($ENV{GENOME_DS_DGIDB_TYPE}) {
    Carp::copak('Environment variable GENOME_DS_DGIDB_TYPE must be set in your environment or by a site configuration module');
}

=pod

=head1 NAME

GENOME_DS_DGIDB_TYPE

=head1 DESCRIPTION

The GENOME_DS_DGIDB_TYPE environment variable indicates the parent class
for the Dgidb data source.  It should be something that inherits from
UR::DataSource::RDBMS.

=cut

1;
