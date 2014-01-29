package Genome::Env::GENOME_USER_EMAIL;
use base 'Genome::Env::Required';

=pod

=head1 NAME

GENOME_USER_EMAIL

=head1 DESCRIPTION

The GENOME_USER_EMAIL environment variable holds the name of the email
address that the LSF job reports are sent to. If the variable is unset then
the email address is constructed using GENOME_EMAIL_DOMAIN and the username.

=cut

1;
