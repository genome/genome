use strict;
use warnings;

use above "Genome";
use Genome::Utility::Test qw(compare_ok abort);
use Test::More tests => 8;

my $class = 'Genome::Model::Tools::Annotate::AddEvsMaf';
use_ok($class);

eval {
    my $data_dir = Genome::Utility::Test->data_dir($class);
    $data_dir = "$data_dir/v1";
    ok(-d $data_dir, "data_dir exists: $data_dir") or abort;

    # check inputs
    my $anno_file = "$data_dir/test.input";
    ok(-s $anno_file, 'input exists: anno_file') or abort;
    my $vcf_file = "$data_dir/ESP6500SI-V2-SSA137.updatedRsIds.snps_indels.vcf";
    ok(-s $vcf_file, 'input exists: vcf_file') or abort;

    # create and execute command
    my $tmp_file = Genome::Sys->create_temp_file_path();
    my $cmd = $class->create(
        anno_file   => $anno_file,
        vcf_file    => $vcf_file,
        output_file => $tmp_file,
    );
    ok($cmd, 'created command') or abort;
    ok($cmd->execute, 'executed command') or abort;

    # check outputs
    ok(-s $tmp_file, 'output_file has size');
    compare_ok($tmp_file, "$data_dir/test.output", 'output_file matched expected');
};
