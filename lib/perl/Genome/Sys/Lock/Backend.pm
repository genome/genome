package Genome::Sys::Lock::Backend;

use strict;
use warnings;

use Mouse::Role;

has 'is_mandatory' => (is => 'ro', isa => 'Bool', required => 1);

requires qw(
    clear_state
    lock
    release_all
    translate_lock_args
    translate_unlock_args
    unlock
);

1;
