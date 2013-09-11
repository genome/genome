use strict;
use warnings;

use above 'Genome';
use Test::More;
use File::Spec;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

my $class = "Genome::Model::Tools::Vcf::CreateCrossSampleVcf::DetermineVarscanOptions";
use_ok($class) || die;

my $cmd = $class->create(
    output_directory => 'some_output_directory',
);

ok($cmd->execute(), 'Successfully executed command');

is(File::Spec->join('some_output_directory', 'varscan_consensus.vcf'),
    $cmd->output_file, "output_file is as expected");
is(1, $cmd->output_vcf, "output_vcf is as expected");

done_testing();
