package Genome::Env::GENOME_TEST_USEDB;

=pod

=head1 NAME

GENOME_TEST_USEDB

=head1 DESCRIPTION

If the GENOME_TEST_USEDB environment variable has a true value, the primary
data source (Genome::DataSource::GMSchema) is reconfigured to use a test
SQLite database next to the currently running program.

=cut

1;
