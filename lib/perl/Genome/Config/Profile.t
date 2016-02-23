#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::Exception;
use Test::More;

use above 'Genome';

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

my $class = 'Genome::Config::Profile';
use_ok($class);

my $ap = Genome::Config::AnalysisProject->create(name => 'Test Project');
_create_ap_infastructure($ap);
my $profile = Genome::Config::Profile->create_from_analysis_project($ap);
for my $rule_model_map ($profile->config_rule_maps){
    ok($rule_model_map->config->status ne 'disabled', 'Config::Profile::Item is not disabled')
}

subtest 'create_from_config_profile_item' => sub{
    plan tests => 6;

    my ($item) = Genome::Config::AnalysisMenu::Item->get(id => -1234)->config_items;
    my $item11 = Genome::Config::Profile::Item->create(
        analysis_project => $ap,
        status => 'disabled',
    );
    ok($item, 'got config item');

    is($item->status('disabled'), 'disabled', 'config item status is disabled');
    throws_ok(
        sub{ $class->create_from_config_profile_item($item); },
        qr/Cannot use a non current config item to create profile/,
        'failed to create profile with disabled config',
    );
    $item->status('inactive');
    $ap->status('Completed');
    throws_ok(
        sub{ $class->create_from_config_profile_item($item); },
        qr/Cannot use a non current config item to create profile/,
        'failed to create profile with non current AnP',
    );
    $ap->status('Pending');

    my $profile;
    lives_ok(sub{ $profile = $class->create_from_config_profile_item($item); }, 'create profile from item');
    isa_ok($profile, $class);

};

done_testing();

sub _create_ap_infastructure {
    my $ap = shift;
    my $menu_item_contents = _get_menu_item_contents();
    my $menu_item_file = _create_file_with_contents($menu_item_contents);
    my $menu_item = Genome::Config::AnalysisMenu::Item->create(
        id => -1234,
        file_path => $menu_item_file,
        name => 'test menu item',
    );

    my $item1 = Genome::Config::Profile::Item->create(
        analysis_menu_item => $menu_item,
        analysis_project => $ap,
        status => 'active',
    );

    my $item3 = Genome::Config::Profile::Item->create(
        analysis_menu_item => $menu_item,
        analysis_project => $ap,
        status => 'active',
    );

    my $item2 = Genome::Config::Profile::Item->create(
        analysis_menu_item => $menu_item,
        analysis_project => $ap,
        status => 'disabled',
    );

    return;
}

sub _create_file_with_contents {
    my $contents = shift;
    my $file_path = join('', Genome::Sys->create_temp_file_path(), '.yml');

    Genome::Sys->write_file($file_path, $contents);

    return $file_path;
}

sub _get_menu_item_contents {
return <<YML
rules:
  sequencing_platform: solexa
  species_name: human
  'sample->extraction_type': genomic dna
  'library->is_bisulfite_converted': 0

models:
  'Genome::Model::ReferenceAlignment':
    - processing_profile_id: 2635769
      annotation_reference_build_id: 124434505
      dbsnp_build_id: 127786607
      reference_sequence_build_id: 106942997
      roi_track_name: 'tiled_region'
      instrument_data_properties:
        region_of_interest_set_name: [target_region_set_name]
        target_region_set_name: [target_region_set_name]
        subject: sample
YML
}
