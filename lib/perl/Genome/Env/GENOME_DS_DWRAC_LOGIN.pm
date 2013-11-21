package Genome::Env::GENOME_DS_DWRAC_LOGIN;

use Carp;

unless ($ENV{GENOME_DS_DWRAC_LOGIN}) {
    Carp::copak('Environment variable GENOME_DS_DWRAC_LOGIN must be set in your environment or by a site configuration module');
}

=pod

=head1 NAME

GENOME_DS_DWRAC_LOGIN

=head1 DESCRIPTION

The GENOME_DS_DWRAC_LOGIN environment variable holds the login name used
to connect to the Dwrac data source.

=cut

1;
