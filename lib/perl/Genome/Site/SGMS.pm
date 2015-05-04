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
        parent_dir => Genome::Config::get('host_lock_dir')));

Genome::Sys::Lock->add_backend('site',
    Genome::Sys::Lock::FileBackend->new(is_mandatory => 1,
        parent_dir => Genome::Config::get('site_lock_dir')));

1;
