#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

my $pkg = 'Genome::InstrumentData::AlignmentResult::Bwamem';
use_ok($pkg);

subtest "parse params" => sub {
    my @invalid_strings = (
        '-o -m -g', 'banana',
        '-c should_be_int',
        '-t should_be_int',
        '-r should_be_float',
        '-c should_be_int',
        '-D should_be_float',
        '-W should_be int',
        '-m should_be_int',
        '-S should_be_null',
        '-P should_be_null',
        '-e should_be_null',
        );

    for my $invalid (@invalid_strings) {
        eval {
            $pkg->_param_string_to_hash($invalid);
        };

        ok($@, "Param string $invalid is invalid");
    }

    my $params = $pkg->_param_string_to_hash("-t 4 -M");
    ok(exists $params->{M}, "-M flag exists");
    ok(exists $params->{t}, "-t flag exists");
    ok(!$params->{M}, "-M doesn't have a value");
    is($params->{t}, 4, "parsed 4 threads");
};

subtest "required rusage" => sub {
    for my $i (4..8) {
        my $rusage = $pkg->required_rusage(
            aligner_params => "-t $i"
            );
        ok($rusage =~ /cpus >= $i .* -n $i/, "threads set correctly (-t $i)");
    }
};

subtest "align" => sub {
    my $ar = make_alignment_result("-t 2");
    ok($ar, "Created alignment result");
};

done_testing();
