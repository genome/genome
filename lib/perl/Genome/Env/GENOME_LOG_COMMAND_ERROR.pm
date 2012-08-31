package Genome::Env::GENOME_LOG_COMMAND_ERROR;
our $default_value = 'default';

=pod

=head1 NAME

GENOME_LOG_COMMAND_ERROR

=head1 DESCRIPTION

Log various failed command metrics when available, like:
message
package
file
subroutine
line
inferred_message
inferred_file
inferred_line
build_id

Set to 1 to log all errors
Set to 0 to log zero errors

=head1 DEFAULT VALUE

'default'

Only log errors when a build id is present and apipe-builder is the user

=cut

1;
