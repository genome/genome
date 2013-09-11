use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Model::TestHelpers qw(create_test_sample);

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

my $class = "Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateOutputDirectory";
use_ok($class) || die;

my $test_sample_name = 'Some Test Sample Name';
create_test_sample($test_sample_name);

my $base_directory = Genome::Sys->create_temp_directory();
my $cmd = $class->create(
    sample_name => $test_sample_name,
    base_directory => $base_directory,
);
ok($cmd->execute(), 'Successfully executed command');
ok(-d $cmd->output_directory, "Found output directory");

done_testing();
