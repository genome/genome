
package Genome::InlineConfig;
use strict;
use warnings;
use Config;

our $DIRECTORY;
sub DIRECTORY {
    unless(defined $DIRECTORY) {
        $DIRECTORY = $INC{"Genome/InlineConfig.pm"};
        $DIRECTORY =~ s/\.pm(\/|)//;
        $DIRECTORY .= ($Config{use64bitall}) ? "64" : "32";
        unless (-d $DIRECTORY) {
            unless(mkdir $DIRECTORY) {
                die "failed to create directory $DIRECTORY: $!";
            }
        }
    }
    return $DIRECTORY;
}

our $CCFLAGS;
sub CCFLAGS {
    unless (defined($CCFLAGS)) {
        $CCFLAGS = `uname -m` =~ /ia64/?'-D_FILE_OFFSET_BITS=64 -m32':'-D_FILE_OFFSET_BITS=64';
    }
    return $CCFLAGS;
}

1;

