package Genome::Env::GENOME_SYS_LOG_LEVEL;

=pod

=head1 NAME

GENOME_SYS_LOG_LEVEL

=head1 DESCRIPTION

The GENOME_SYS_LOG_LEVEL environment variable can be set to any of the following which are not n/a:

  0  debug       debug_message (and those below) will go to syslog
  1  info        status_message (and those below) will go to syslog
  2  notice      n/a  
  3  warning     warning_message (and those below) will go to syslog
  4  error       error_message (and those below) will got to syslog
  5  critical    n/a 
  6  alert       n/a
  7  emergency   n/a
  

=head1 DEFAULT VALUE

No value is set by default, meaning no syslog entries will be made.

=cut

1;

