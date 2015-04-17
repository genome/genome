package Genome::Site::SGMS;
use strict;
use warnings;

BEGIN {
    unless ($INC{'Genome/Sys/Lock.pm'}) {
        use Genome::Sys::Lock;
    }
};

require Genome::Sys::Lock::FileBackend;

Genome::Sys::Lock->add_backend('host',
    Genome::Sys::Lock::FileBackend->new(is_mandatory => 1,
        parent_dir => $ENV{GENOME_HOST_LOCK_DIR}));

Genome::Sys::Lock->add_backend('site',
    Genome::Sys::Lock::FileBackend->new(is_mandatory => 1,
        parent_dir => $ENV{GENOME_SITE_LOCK_DIR}));

1;
