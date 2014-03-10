
#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Test::More;

my $class = 'Genome::Config::AnalysisProject::Command::AddConfigFile';
use_ok($class);


my $file = Genome::Sys->create_temp_file_path();
my $contents = 'Test Contents';
Genome::Sys->write_file($file, $contents);


my $analysis_project = Genome::Config::AnalysisProject->create(
    name => 'Test Project'
);

my $cmd = $class->create(
    analysis_project => $analysis_project,
    config_file => $file,
);

$cmd->execute();
my $config_profile_item = Genome::Config::Profile::Item->get(analysis_project => $analysis_project);

ok($config_profile_item, 'it should create a config profile item and associate it with the analysis project');
ok($config_profile_item->_is_concrete, 'it should be concrete');
is(Genome::Sys->read_file($config_profile_item->file_path), $contents, 'it should copy the file to the allocation of the profile item');

done_testing();

1;
