package Genome::Site::SGMS;
use strict;
use warnings;
use File::Spec;

BEGIN {
    unless ($INC{'Genome/Sys/Lock.pm'}) {
        use Genome::Sys::Lock;
    }
};

require Genome::Sys::Lock::FileBackend;
Genome::Sys::Lock->add_backend('site',
    Genome::Sys::Lock::FileBackend->new(is_mandatory => 1,
        parent_dir => File::Spec->catdir(
            $ENV{GENOME_LOCK_DIR}, 'site')));

Genome::Sys::Lock->add_backend('unknown',
    Genome::Sys::Lock::FileBackend->new(is_mandatory => 1,
        parent_dir => File::Spec->catdir(
            $ENV{GENOME_LOCK_DIR}, 'unknown')));

Genome::Sys::Lock->add_backend('tgisan',
    Genome::Sys::Lock::FileBackend->new(is_mandatory => 1,
        parent_dir => File::Spec->catdir(
            $ENV{GENOME_LOCK_DIR}, 'tgisan')));
1;
