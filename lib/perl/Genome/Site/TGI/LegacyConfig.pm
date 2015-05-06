package Genome::Config;

use strict;
use warnings;

our $VERSION = $Genome::VERSION;

use UR;
use Sys::Hostname;

my $arch_os;
sub arch_os {
    unless ($arch_os) {
        $arch_os = `uname -m`;
        chomp($arch_os);
    }
    return $arch_os;
}

# in dev mode we use dev search, dev wiki, dev memcache, etc, but production database still ;)
my $dev_mode = ( Genome::Config::get('dev_mode') || UR::DBI->no_commit );
if ($dev_mode) {
    my $h = hostname;
    warn "***** GENOME_DEV_MODE ($h) *****";
}

sub dev_mode {
    shift;
    if (@_ && !Genome::Config::get('dev_mode')) {
        $dev_mode = shift;
    }

    return $dev_mode;
}

sub auth_user {

    my ($class, $u) = @_;
    my $auth_user = Genome::Sys->username();
    if (defined($u)) {
        $auth_user = $u;
    }

    return $auth_user;
}

sub reference_sequence_directory {
    return join('/', Genome::Config::get('model_root'), 'reference_sequences');
}

1;

=pod

=head1 NAME

Genome::Config - environmental configuration for the genome modeling tools

=head1 DESCRIPTION

The methods in this module are undergoing heavy refactoring and should be ignored until a later release.

=head1 AUTHORS

This software is developed by the analysis and engineering teams at 
The Genome Center at Washington Univiersity in St. Louis, with funding from 
the National Human Genome Research Institute.

=head1 LICENSE

This software is copyright Washington University in St. Louis.  It is released under
the Lesser GNU Public License (LGPL) version 3.  See the associated LICENSE file in
this distribution.

=head1 BUGS

For defects with any software in the genome namespace,
contact genome-dev@genome.wustl.edu.

=cut

