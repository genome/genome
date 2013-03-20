use strict;
use warnings;

use above "Genome";
use Genome::Utility::Test qw(compare_ok);
use Test::More tests => 8;

my $class = 'Genome::Model::Tools::Annotate::AddRsid';
use_ok($class);

do {
    my $data_dir = Genome::Utility::Test->data_dir($class);
    ok(-d $data_dir, "data_dir exists: $data_dir") or die;

    # check inputs
    my $anno_file = "$data_dir/test.annotate.top2";
    ok(-s $anno_file, 'input exists: anno_file') or die;
    my $vcf_file = "$data_dir/test.vcf.gz";
    ok(-s $vcf_file, 'input exists: vcf_file') or die;

    # create and execute command
    my $tmp_file = Genome::Sys->create_temp_file_path();
    my $cmd = $class->create(
        anno_file   => $anno_file,
        vcf_file    => $vcf_file,
        output_file => $tmp_file,
    );
    ok($cmd, 'created command') or die;
    ok($cmd->execute, 'executed command') or die;

    # check outputs
    ok(-s $tmp_file, 'output_file has size');
    compare_ok($tmp_file, "$data_dir/test.out", 'output_file matched expected');
};
