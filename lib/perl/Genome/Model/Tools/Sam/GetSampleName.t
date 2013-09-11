#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Model::Tools::Dindel::TestHelpers qw(
    get_test_dir
);

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

my $class = 'Genome::Model::Tools::Sam::GetSampleName';
use_ok($class);

my $VERSION = 1; # Bump this each time test data changes

my $test_dir = get_test_dir($class, $VERSION);
my $input_bam = File::Spec->join($test_dir, '134053361-region-limited.bam');
ok(-s $input_bam, "Found input Bam file");

my $cmd = $class->create(bam_file => $input_bam);
ok($cmd->execute(), "Successfully executed command");
is('H_LH-ED0001-58404353', $cmd->sample_name, 'Found the correct sample name');

done_testing();
