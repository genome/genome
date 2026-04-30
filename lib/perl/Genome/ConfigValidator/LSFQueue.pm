package Genome::ConfigValidator::LSFQueue;

use strict;
use warnings;

use Mouse;

with qw(Genome::ConfigValidatorBase);

sub check {
    my ($self, $value) = @_;
    my $jdb = Genome::Config::get('job_dispatch_backend');
    if ($jdb eq 'lsf') {
        system(qq(bqueues $value 1> /dev/null 2>&1));
        return $? == 0;
    } elsif ($jdb eq 'slurm') {
        system(qq(scontrol show partition $value 1> /dev/null 2>&1));
        return $? == 0;
    } else {
        die 'unknown job dispatch backend';
    }
}

sub message {
    return 'an LSF queue';
}

1;
