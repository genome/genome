package Genome::Sys::Lock;

use strict;
use warnings;

BEGIN {
    unless ($INC{'Genome/Sys/Lock.pm'}) {
        die 'must load Genome::Sys::Lock first';
    }
};


require Genome::Sys::Lock::FileBackend;
Genome::Sys::Lock->add_backend('site',
    Genome::Sys::Lock::FileBackend->new(is_mandatory => 1,
        parent_dir => $ENV{GENOME_LOCK_DIR}));
Genome::Sys::Lock->add_backend('tgisan',
    Genome::Sys::Lock::FileBackend->new(is_mandatory => 1,
        parent_dir => $ENV{GENOME_FS_LOCAL_NETWORK_CACHE}));

Genome::Sys::Lock->add_backend('unknown',
    Genome::Sys::Lock::FileBackend->new(is_mandatory => 1,
        parent_dir => '/'));


if ($ENV{GENOME_NESSY_SERVER}) {
    require Genome::Sys::Lock::NessyBackend;
    my $is_mandatory = $ENV{GENOME_NESSY_MANDATORY} ? 1 : 0;
    my $nessylock = Genome::Sys::Lock::NessyBackend->new(
        url => $ENV{GENOME_NESSY_SERVER},
        is_mandatory => $is_mandatory,
    );
    Genome::Sys::Lock->add_backend('site', $nessylock);

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
