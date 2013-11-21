package Genome::Env::GENOME_DS_GMSCHEMA_SERVER;

use Carp;

unless ($ENV{GENOME_DS_GMSCHEMA_SERVER}) {
    Carp::copak('Environment variable GENOME_DS_GMSCHEMA_SERVER must be set in your environment or by a site configuration module');
}

=pod

=head1 NAME

GENOME_DS_GMSCHEMA_SERVER

=head1 DESCRIPTION

The GENOME_DS_GMSCHEMA_SERVER environment variable holds database server
connection details for the GMSchema database.  Its value is used to build
the DBI connection string after the DBI driver name.

=cut

1;
