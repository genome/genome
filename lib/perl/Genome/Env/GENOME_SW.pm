package Genome::Env::GENOME_SW;
our $default_value = '/var/lib/genome/sw';

=pod

=head1 NAME

GENOME_SW

=head1 DESCRIPTION

The GENOME_SW environment variable is the root directory for data base storage.

=head1 DEFAULT VALUE

The Genome::Sys module typically leans on OS package management where advanced package management exists (Debian/Ubuntu, RedHat, etc.)

This path is a fallback for systems with locally produced software which is not packaged:

 /var/lib/genome/sw

=cut

1;
