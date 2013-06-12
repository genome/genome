# this module patches some parts of the source tree to it functional for release 
package Genome::Patch;
use strict;
use warnings;
use Sub::Install;

# the formula for handling Illumina fastqs is wrong and needs to be fixed
# on the master branch, compensate for now:
require Genome::InstrumentData::Solexa;
my $rqc_old = \&Genome::InstrumentData::Solexa::resolve_quality_converter;
my $rqc_new = \&Genome::InstrumentData::Solexa::NEW_resolve_quality_converter;
delete $Genome::InstrumentData::Solexa::{resolve_quality_converter};
delete $Genome::InstrumentData::Solexa::{NEW_resolve_quality_converter};
Sub::Install::install_sub({
    code => $rqc_new, 
    into => "Genome::InstrumentData::Solexa",
    as   => 'resolve_quality_converter',
});

# Errors in the base class are redundant with db constraints
# and cause issues because some data types need to be updated
# with the pg transition.  Supress for now.
require UR;
my $orig_errors = \&UR::Object::__errors__;
sub __patched_errors__ { return(); };
delete $UR::Object::{__errors__};
Sub::Install::install_sub({
    code => \&__patched_errors__,
    into => "UR::Object",
    as   => '__errors__',
});

# the alignment summary executable is an ancient C++ hack not compiled for release yet
require Genome::InstrumentData::AlignmentResult;
delete $Genome::InstrumentData::AlignmentResult::{_use_alignment_summary_cpp};
Sub::Install::install_sub({
    code => sub { 0 },
    into => 'Genome::InstrumentData::AlignmentResult',
    as   => '_use_alignment_summary_cpp',
});

1;

