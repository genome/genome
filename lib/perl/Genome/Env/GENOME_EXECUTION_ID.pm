package Genome::Env::GENOME_EXECUTION_ID;

use UR;
sub default_value { UR::Object::Type->autogenerate_new_object_id_uuid() }

=pod

=head1 NAME

GENOME_EXECUTION_ID

=head1 DESCRIPTION

A UUID value representing one case of executing a Genome program.

=cut

1;
