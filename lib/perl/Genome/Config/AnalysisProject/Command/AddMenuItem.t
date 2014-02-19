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

my $inst_data_1 = Genome::InstrumentData::Imported->create();
my $inst_data_2 = Genome::InstrumentData::Imported->create();

my $ap = Genome::Config::AnalysisProject->create(
    name => 'Test Project'
);

my $b1 = Genome::Config::AnalysisProject::InstrumentDataBridge->create( analysis_project => $ap, instrument_data => $inst_data_1, status => 'processed');

my $b2 = Genome::Config::AnalysisProject::InstrumentDataBridge->create( analysis_project => $ap, instrument_data => $inst_data_2, status => 'processed');

my $item = Genome::Config::AnalysisMenu::Item->create( name => 'test_item', file_path => '/tmp/idontexist');

my $cmd = $class->create(
    analysis_project => $ap,
    analysis_menu_items => [$item],
    reprocess_existing => 1,
);

$cmd->execute();

ok(Genome::Config::Profile::Item->get(analysis_project => $ap, analysis_menu_item => $item), 'added item to the project');

ok($b1->status eq 'new', 'first instrument data set to new');
ok($b2->status eq 'new', 'second instrument data set to new');

done_testing();
