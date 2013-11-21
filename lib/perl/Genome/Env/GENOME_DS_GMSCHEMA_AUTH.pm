package Genome::Env::GENOME_DS_GMSCHEMA_AUTH;

use Carp;

unless ($ENV{GENOME_DS_GMSCHEMA_AUTH}) {
    Carp::copak('Environment variable GENOME_DS_GMSCHEMA_AUTH must be set in your environment or by a site configuration module');
}

=pod

=head1 NAME

GENOME_DS_GMSCHEMA_AUTH

=head1 DESCRIPTION

The GENOME_DS_GMSCHEMA_AUTH environment variable holds the login password used
to connect to the GMSchema data source.

=cut

1;
