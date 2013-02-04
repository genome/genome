
package Genome::UrPatch;

# this module monkey patches UR
# to get around some things in the releaseable GMS.

use UR; 
use Sub::Install;

my $orig_errors = \&UR::Object::__errors__;

sub __patched_errors__ {
    $orig_errors->(@_);
    # Errors in the base class are redundant with db constraints
    # and cause issues because some data types need to be updated
    # with the pg transition.  Supress for now.
    return();
};

no warnings;
Sub::Install::install_sub({
    code => \&Genome::UrPatch::__patched_errors__,
    into => "UR::Object",
    as   => '__errors__',
});

1;

