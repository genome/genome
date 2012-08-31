#!/usr/bin/env genome-perl
use strict;
use warnings;
use File::Basename;
chdir File::Basename::dirname(__FILE__);
my @dirs = (glob("*/debian/*"), glob("*/ubuntu-lucid/*"), glob("vendor/*/debian/*"), glob("vendor/*/ubuntu-lucid/*"));

@dirs = grep { 
        if (m{^(vendor/|)(.*)/(debian|ubuntu-lucid)/(.*)}) {
            if ($2 eq $4) {
                ($_);
            }
            else {
                #print "missed $_\n";
                ();
            }
        }
        else {
            die "odd: $_";
        }
    }
    @dirs;

if (@dirs == 0) {
    print "No directories with debian build cruft found under $ENV{PWD} (checkted for 'debian' and 'ubuntu-lucid').\n";
    exit;
}

print "Removing all data from these directories under $ENV{PWD}:\n\t",join("\n\t",@dirs),"\n\n";
print "Press ctr-c in 5 seconds if you don't want to do this...\n";

sleep 5;

for my $dir (@dirs) {
    my $cmd = "rm -rf '$dir'";
    print "RUN: $cmd\n";
    my $rv = system $cmd;
    $rv /= 256;
    if ($rv) {
        warn "error for $dir: $!";
        next;
    }
}

