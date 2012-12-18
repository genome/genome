#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

# Use
use_ok('Genome::Model::Tools::Sx::Functions') or die;

# Quality calculations
is(Genome::Model::Tools::Sx::Functions->calculate_quality('BBAB<BBBBAB=??#@?8@1(;>A::(4@?--98#########################################'), 960, 'calculate quality'); 
is(Genome::Model::Tools::Sx::Functions->calculate_average_quality('BBAB<BBBBAB=??#@?8@1(;>A::(4@?--98#########################################'), 13, 'calculate average quality'); 
is(Genome::Model::Tools::Sx::Functions->calculate_qualities_over_minumum('BBAB<BBBBAB=??#@?8@1(;>A::(4@?--98#########################################', 20), 27, 'calculate qualities over min'); 
is(Genome::Model::Tools::Sx::Functions->minimum_quality('BBBBBBBBBB<BBBBBBBBBB'), 27, 'minimum quality'); 
is(Genome::Model::Tools::Sx::Functions->maximum_quality('BBBBBBBBBBCBBBBBBBBBB'), 34, 'maximum quality'); 

# Config processing
# fail
ok(!Genome::Model::Tools::Sx::Functions->hash_to_config, 'Hash to config fails w/o hash');
is(Genome::Model::Tools::Sx::Functions->error_message, 'No hash to convert to config!', 'Correct error message');
ok(!Genome::Model::Tools::Sx::Functions->hash_to_config(array => []), 'Hash to config fails w/ array value');
like(Genome::Model::Tools::Sx::Functions->error_message, qr/^Unsupported data type \(ARRAY\) in hash!/, 'Correct error message');

ok(!Genome::Model::Tools::Sx::Functions->config_to_hash, 'Config to hash fails w/o hash');
is(Genome::Model::Tools::Sx::Functions->error_message, 'No config to convert to hash!', 'Correct error message');
ok(!Genome::Model::Tools::Sx::Functions->config_to_hash('array=val1:array=val2'), 'Config to hash fails w/ multiple values');
is(Genome::Model::Tools::Sx::Functions->error_message, 'Duplicate key (array) in config! array=val1:array=val2', 'Correct error message');

# success
my %expected_hash = (
    name => 'Barack',
    type => 'president',
    other => 'Obama'
);
my $expected_config = 'name=Barack:other=Obama:type=president';
my $config = Genome::Model::Tools::Sx::Functions->hash_to_config(%expected_hash);
is($config, $expected_config, 'Hash to config');
my %hash = Genome::Model::Tools::Sx::Functions->config_to_hash($config);
is_deeply(\%hash, \%expected_hash, 'Config to hash');

done_testing();
