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

is($cmd->output_file,
    File::Spec->join('some_output_directory', 'varscan_consensus.vcf.gz'),
    "output_file is as expected");
is($cmd->output_vcf, 1, "output_vcf is as expected");

done_testing();
