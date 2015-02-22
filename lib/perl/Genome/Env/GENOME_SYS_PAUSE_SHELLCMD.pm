package Genome::Env::GENOME_SYS_PAUSE_SHELLCMD;
sub default_value { '/var/lib/genome/testsuite-inputs' }

=pod

=head1 NAME

GENOME_SYS_PAUSE_SHELLCMD

=head1 DESCRIPTION

When a shell command is to be executed with a specific pattern, creates a temp file instead with that command, and waits for it to be deleted.

=head1 DEFAULT VALUE

 not set

=cut

1;

