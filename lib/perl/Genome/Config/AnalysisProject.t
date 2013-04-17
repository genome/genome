#!/usr/bin/env genome-perl
use strict;
use warnings;

$ENV{UR_DBI_NO_COMMIT} = 1;

use Test::More;
use above "Genome";
use Carp::Always;

my $class = 'Genome::Config::AnalysisProject';

use_ok($class);

my $config = Genome::Config::Set->create();
my $analysis_config_set = Genome::Config::Set->create();

my $menu_item = Genome::Config::AnalysisMenuItem->create(
    name => 'Test menu item',
    configuration_set => $analysis_config_set,
);

my $test_obj = Genome::Config::AnalysisProject->create(
    created_by => Genome::Sys->username,
    name => 'Test Project',
    _analysis_menu_item => $menu_item,
    _configuration_set => $config,
);

my $config_reader = $test_obj->get_configuration_reader;
ok($config_reader, 'it constructs and returns the reader object when called');
isa_ok($config_reader, 'Genome::Config::MaskedConfigurationReader', 'it returns the appropriate object type');
ok($test_obj->configuration_reader, 'it memoizes the configuration reader object');

done_testing();


