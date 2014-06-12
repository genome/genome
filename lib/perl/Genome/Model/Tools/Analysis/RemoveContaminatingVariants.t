use strict;
use warnings;

use above "Genome";
use Genome::Utility::Test qw(compare_ok abort);
use Test::More tests => 7;

my $class = 'Genome::Model::Tools::Analysis::RemoveContaminatingVariants';
use_ok($class);

eval {
    my $data_dir = Genome::Utility::Test->data_dir($class);
    $data_dir = "$data_dir/v1";
    print STDER $data_dir . "\n";
    ok(-d $data_dir, "data_dir exists: $data_dir") or abort;

    # check inputs
    my $infile = "$data_dir/testdata";
    ok(-s $infile, 'input exists') or abort;
    
    # create and execute command
    my $tmp_file = Genome::Sys->create_temp_file_path();
    my $cmd = $class->create(
        input_file   => $infile,
        output_file => $tmp_file,
    );
    ok($cmd, 'created command') or abort;
    ok($cmd->execute, 'executed command') or abort;

    # check outputs
    ok(-s $tmp_file, 'output_file has size');
    compare_ok($tmp_file, "$data_dir/testdata.out", 'output_file matched expected');
};
