package Genome::Env::GENOME_DS_OLTP_TYPE;
use base 'Genome::Env::Required';

=pod

=head1 NAME

GENOME_DS_OLTP_TYPE

=head1 DESCRIPTION

The GENOME_DS_OLTP_TYPE environment variable indicates the parent class
for the Oltp data source.  It should be something that inherits from
UR::DataSource::RDBMS.

=cut

1;
