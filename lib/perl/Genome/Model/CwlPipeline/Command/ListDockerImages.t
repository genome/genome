#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Genome::Test::Factory::ProcessingProfile::CwlPipeline;

use Test::More tests => 5;

my $class = 'Genome::Model::CwlPipeline::Command::ListDockerImages';

use_ok($class);

my $pp = Genome::Test::Factory::ProcessingProfile::CwlPipeline->setup_object(
    main_workflow_file => 'https://raw.githubusercontent.com/genome/analysis-workflows/eb0092603bf57acb7bda08a06e4f2f1e2a8c9b6d/definitions/pipelines/somatic_exome.cwl',
);

my $output = Genome::Sys->create_temp_file_path;

my $cmd = $class->create(
    processing_profile => $pp,
    output_file => $output,
);
isa_ok($cmd, $class, 'created command');

SKIP: {
    skip "needs cwltool", 3 unless -x '/usr/local/bin/cwltool';
    ok($cmd->execute, 'executed command');

    my @data = Genome::Sys->read_file($output);

    ok($data[0] !~ /^ /, 'first line does not start with image');

    my @image_lines = grep { $_ =~ /^  / } @data;
    ok(scalar(@image_lines) > 0, 'theoretically included some docker images');
}

