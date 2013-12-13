package Genome::Env::GENOME_DS_GMSCHEMA_SERVER;
use base 'Genome::Env::Required';

=pod

=head1 NAME

GENOME_DS_GMSCHEMA_SERVER

=head1 DESCRIPTION

The GENOME_DS_GMSCHEMA_SERVER environment variable holds database server
connection details for the GMSchema database.  Its value is used to build
the DBI connection string after the DBI driver name.

=cut

1;
