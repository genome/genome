#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More tests => 4;

use above 'Genome';

use Genome::Test::Factory::AnalysisProject;

my $class = 'Genome::Config::AnalysisProject::Command::PossibleVolumes';

use_ok($class);

my $dg = Genome::Disk::Group->__define__(
    name => 'testing_possible_volumes',
);
my $dv1 = Genome::Disk::Volume->__define__(
    mount_path => '/tmp/testing_possible_volumes',
    can_allocate => 1,
);

my $da1 = Genome::Disk::Assignment->__define__(
    group_id => $dg->id,
    volume_id => $dv1->id,
);

my $dv2 = Genome::Disk::Volume->__define__(
    mount_path => '/tmp/testing_possible_volumes2',
    can_allocate => 1,
);

my $da2 = Genome::Disk::Assignment->__define__(
    group_id => $dg->id,
    volume_id => $dv2->id,
);
my $anp = Genome::Test::Factory::AnalysisProject->setup_object;

my $cmd = $class->create(analysis_project => $anp);

my $env_file = Genome::Sys->create_temp_file_path;
Genome::Sys->write_file($env_file, <<EOENV
disk_group_alignments: testing_possible_volumes
disk_group_models: testing_possible_volumes
disk_group_scratch: testing_possible_volumes
EOENV
);

Genome::Config::AnalysisProject::Command::AddEnvironmentFile->execute(
    environment_file => $env_file,
    analysis_project => $anp
);



isa_ok($cmd, $class, 'created command');

my ($fh, $file) = Genome::Sys->create_temp_file();
my $rv;
{
    local *STDOUT = $fh;
    $rv = $cmd->execute;
    close *STDOUT;
}

ok($rv, "executed command");
my $expected = <<EOEXPECTED
/tmp/testing_possible_volumes
/tmp/testing_possible_volumes2
EOEXPECTED
;

my $diff = Genome::Sys->diff_file_vs_text($file, $expected);
ok(!$diff, 'output matched expected');
