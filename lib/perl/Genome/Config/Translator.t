#!/usr/bin/env genome-perl
use strict;
use warnings;

$ENV{UR_DBI_NO_COMMIT} = 1;

use Test::More;
use above "Genome";

my $class = 'Genome::Config::Translator';
use_ok($class);

eval_result(Genome::Config::Translator::_wrap_single_model_hashes(_single_model_type()));
eval_result(Genome::Config::Translator::_wrap_single_model_hashes(_multiple_model_types()));

done_testing();

sub eval_result {
    my $result = shift;
    for (values %$result) {
        is('ARRAY', ref $_, 'Expect all results to be array refs');
    }
}

sub _single_model_type {
    return {
        ModelType1 => _one_config()
    };
}

sub _multiple_model_types {
    return {
        ModelType1 => _one_config(),
        ModelType2 => _multiple_configs(),
    };
}

sub _one_config {
    return {
        turkey => 'sandwich',
    };
}

sub _multiple_configs {
    return [
        _one_config(),
        _one_config(),
    ];
}