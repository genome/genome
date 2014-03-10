package Genome::Env::GENOME_EMAIL_SMTP_SERVER;
use base 'Genome::Env::Required';

=pod

=head1 NAME

GENOME_EMAIL_SMTP_SERVER

=head1 DESCRIPTION

This variable is used to set the SMTP server used to send email.

=head1 DEFAULT VALUE

$(HOSTNAME)
=======
The GENOME_EMAIL_SMTP_SERVER environment variable holds the name of the email
server used to send email through.

=cut

1;
