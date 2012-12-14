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

done_testing();
