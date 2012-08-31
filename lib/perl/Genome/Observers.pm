package Genome::Observers;

use strict;
use warnings;
use File::Basename;

# Observers for Genome classes should go in Genome/Observers. These
# observers exist to separate dependencies between classes. Each module
# should following the naming schema <class to observe><related class>.
# For example, a module that contains observers for the relationship 
# between Samples and Libraries should be called SampleLibrary.pm.

# Try to avoid cascading deletes in these observers unless you are 
# absolutely certain it's the correct thing to do. In most cases, you 
# should throw an exception rather than cascade.

# Loads all observers in Genome/Observers/*
my $base_path = __FILE__;
$base_path =~ s/.pm$//;
my @paths = glob("$base_path/*");

for my $path (@paths) {
    my $file = basename($path);
    $file =~ s/.pm$//;

    my $package = __PACKAGE__ . "::$file";
    eval "use $package";
    die $@ if $@;
}

1;

