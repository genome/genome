#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use Sub::Override;

use Test::More tests => 3;

use above 'Genome';

my $class = 'Genome::Ptero::ExecutionInfo';

use_ok($class);

my $execution = &_fake_execution;
my $job_info = &_fake_job_info;

my $override = Sub::Override->new('Ptero::HTTP::make_request_and_decode_response', sub {
    my %args = @_;
    if ($args{url} eq $execution->{data}{jobUrl}) {
        return $job_info;
    } else {
        die 'Unexpected request to Ptero server!';
    }
});

my $cmd = $class->create(
    execution_id => -1,
    display_inputs => 1,
);
isa_ok($cmd, $class, 'created command');

$cmd->__execution($execution);

ok($cmd->execute, 'executed command');


sub _fake_execution {
    my $e = {
        inputs => { test => 1 },
        data => {
            jobUrl => 'http://lsf.example.com/v1/jobs/12345678-90ab-cdef-1234-567890abcdef',
        },
        status => 'running',
    };

    bless $e, 'Ptero::Proxy::Workflow::Execution';
    return $e;
}

sub _fake_job_info {
    my $info = {
        lsfJobId => '12345678',
        stdout => "Here is some output.\n",
        stderr => "Here is an error.\n",
        status => 'running',
        jobId => '12345678-90ab-cdef-1234-567890abcdef',
        command => 'echo This is fake.',

    };

    return $info;
}
