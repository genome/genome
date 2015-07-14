#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Test::More;
use Genome::Test::Factory::InstrumentData::Solexa;

my $class = 'Genome::Config::AnalysisProject::Command::AddMenuItem';
use_ok($class);

my $inst_data_1 = Genome::Test::Factory::InstrumentData::Solexa->setup_object;
my $inst_data_2 = Genome::Test::Factory::InstrumentData::Solexa->setup_object;

my $ap = Genome::Config::AnalysisProject->create(
    name => 'Test Project',
    run_as => 'nobody',
);

my $b1 = Genome::Config::AnalysisProject::InstrumentDataBridge->create( analysis_project => $ap, instrument_data => $inst_data_1, status => 'processed');

my $b2 = Genome::Config::AnalysisProject::InstrumentDataBridge->create( analysis_project => $ap, instrument_data => $inst_data_2, status => 'processed');

my $item = Genome::Config::AnalysisMenu::Item->create( name => 'test_item', file_path => '/tmp/idontexist', description => 'test');

my $tag = Genome::Config::Tag->create( name => 'test tag for add-menu-item' );

my $cmd = $class->create(
    analysis_project => $ap,
    analysis_menu_items => [$item],
    reprocess_existing => 1,
    tags => [$tag],
);

$cmd->execute();

my $profile_item = Genome::Config::Profile::Item->get(analysis_project => $ap, analysis_menu_item => $item);
ok($profile_item, 'added item to the project');

is($profile_item->tags, $tag, 'tag associated with item');

ok($b1->status eq 'new', 'first instrument data set to new');
ok($b2->status eq 'new', 'second instrument data set to new');

ok(UR::Context->commit(), 'created objects are valid');

$ap->status('Deprecated');
my $fail_cmd = $class->create(
    analysis_project => $ap,
    analysis_menu_items => [$item],
);
ok(!$fail_cmd->execute,'fail command on project with invalid status');

done_testing();
