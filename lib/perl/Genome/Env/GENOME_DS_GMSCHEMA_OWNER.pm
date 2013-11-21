package Genome::Env::GENOME_DS_GMSCHEMA_OWNER;

use Carp;

unless ($ENV{GENOME_DS_GMSCHEMA_OWNER}) {
    Carp::copak('Environment variable GENOME_DS_GMSCHEMA_OWNER must be set in your environment or by a site configuration module');
}

=pod

=head1 NAME

GENOME_DS_GMSCHEMA_OWNER

=head1 DESCRIPTION

The GENOME_DS_GMSCHEMA_OWNER environment variable holds the default
schema/owner for the tables in the GMSchema data source.

=cut

1;
