#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 4;

use_ok('Genome::Model::Tools::Picard::Downsample') or die;

my $cmd = Genome::Model::Tools::Picard::Downsample->create(
    input_file => 'blah',
    output_file => 'blah',
    downsample_ratio => 0,
);
ok($cmd, 'create downsample command w/ invalid downsample_ratio');
my @errors = $cmd->__errors__;
is(@errors, 1, 'correct number of errors');
is($errors[0]->__display_name__, "INVALID: property 'downsample_ratio': Must be greater than 0 and less than 1! 0", 'correct error __display_name__');

done_testing();
