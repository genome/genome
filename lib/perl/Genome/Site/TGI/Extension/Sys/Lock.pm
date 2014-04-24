package Genome::Site::TGI::Extension::Sys::Lock;

use strict;
use warnings;

use Genome::Sys::Lock;

sub import {
    Genome::Sys::Lock->add_backend('Genome::Sys::NessyLock', 0);
}

1;
