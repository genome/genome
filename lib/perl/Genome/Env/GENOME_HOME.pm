package Genome::Env::GENOME_HOME;
use base 'Genome::Env::Required';

=pod

=head1 NAME

GENOME_HOME

=head1 DESCRIPTION

The GENOME_HOME is the directory of the GMS install, typically /opt/gms/$GENOME_SYS_ID.

Other GMS instances may be mounted next-to the home install, in the same directory.

=head1 DEFAULT VALUE

/opt/gms/$GENOME_SYS_ID

=cut

1;
