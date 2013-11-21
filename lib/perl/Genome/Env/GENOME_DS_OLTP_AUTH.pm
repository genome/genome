package Genome::Env::GENOME_DS_OLTP_AUTH;

use Carp;

unless ($ENV{GENOME_DS_OLTP_AUTH}) {
    Carp::copak('Environment variable GENOME_DS_OLTP_AUTH must be set in your environment or by a site configuration module');
}

=pod

=head1 NAME

GENOME_DS_OLTP_AUTH

=head1 DESCRIPTION

The GENOME_DS_OLTP_AUTH environment variable holds the login password used
to connect to the Oltp data source.

=cut

1;
