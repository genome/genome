package Genome::Env::GENOME_HOME;
our $default_value = '/opt/gms/' $ENV{GENOME_SYS_ID};

=pod

=head1 NAME

GENOME_HOME

=head1 DESCRIPTION

The GENOME_HOME environment variable is the root directory for a given GMS instance.
By convention it is under /opt/gms/, in a subdirectory for the giving GMS specified by its GENOME_SYS_ID.

This allows interaction between GMS systems to occur at a simple filesystem level, depending on permissions.

=head1 DEFAULT VALUE

 /opt/gms/$GENOME_SYS_ID/

=cut

1;
