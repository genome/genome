#!/usr/bin/env genome-perl
use strict;
use warnings;

use Test::More;
use above "Genome";
use Genome::Utility::Test;
use Carp::Always;

my $class = 'Genome::Config::MaskedConfigurationReader';

use_ok($class);
my $data_dir = Genome::Utility::Test->data_dir($class, 1);
ok(-d $data_dir, "data_dir exists: $data_dir") or die;

my $tmpdir = Genome::Sys->create_temp_directory();

`cp -r $data_dir/config $tmpdir/`;

my $config_reader = Genome::Config::MaskedConfigurationReader->create(
    config_handler              => Genome::Config::TreeHandler->create(base_path => "$tmpdir/config"),
    mask_handler                => Genome::Config::TreeHandler->create(base_path => "$data_dir/mask"),
    default_handler             => Genome::Config::TreeHandler->create(base_path => "$data_dir/default"),
    configuration_parser        => Genome::Config::Parser::YAML->create(),
    configuration_copy_strategy => Genome::Config::CopyStrategy::TreeCopy->create(),
);

my $result1 = $config_reader->get_config(taxon => 'human', type => 'dna');
ok($result1->{already_here}, 'do not overwrite config if it has already been copied');

my $result2 = $config_reader->get_config(taxon => 'human', type => 'rna');
ok($result2, 'got config that was allowed by the mask, in the default, but not copied yet');

my $result3 = $config_reader->get_config(taxon => 'mouse', type => 'rna');
ok(!$result3, 'dont allow config if it isnt in the mask');

done_testing();
