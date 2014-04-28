package Genome::Sys::Lock;

use strict;
use warnings;

BEGIN {
    unless ($INC{'Genome/Sys/Lock.pm'}) {
        die 'must load Genome::Sys::Lock first';
    }
};

if ($ENV{GENOME_NESSY_SERVER}) {
    require Genome::Sys::Lock::NessyBackend;
    my $nessylock = Genome::Sys::Lock::NessyBackend->new(
        url => 'http://nessy.gsc.wustl.edu/',
        is_mandatory => 0,
    );
    push @Genome::Sys::Lock::backends, $nessylock;

    UR::Context->process->add_observer(
        aspect => 'sync_databases',
        callback => sub {
            my ($ctx, $aspect, $sync_db_result) = @_;
            if ($sync_db_result) {
                use vars '@CARP_NOT';
                local @CARP_NOT = (@CARP_NOT, 'UR::Context');
                foreach my $claim ($nessylock->claims) {
                    $claim->validate
                        || Carp::croak(sprintf('Claim %s failed to verify during commit', $claim->resource_name));
                }
            }
        }
    );
}

1;
