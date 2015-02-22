package Genome::Env::GENOME_TEST_FILLDB;

=pod

=head1 NAME

GENOME_TEST_FILLDB

=head1 DESCRIPTION

If the GENOME_TEST_FILLDB environment variable has a true value, a test
SQLite database is created next to the currently running program, and all
data loaded from the primary database is copied into this SQLite database.

=cut

1;
