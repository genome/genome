#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);
use File::Spec;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{NO_LSF} = 1;
};

my $class = 'Genome::Model::Tools::Dindel::MakeDindelWindows';
use_ok($class);

my $VERSION = 0; # Bump this each time tests data changes

my $test_dir = File::Spec->join(Genome::Utility::Test->data_dir($class), "v$VERSION");
diag "Test data located at $test_dir\n";

my $input_file = File::Spec->join($test_dir, 'input.dindel');
ok(-s $input_file, 'Found input dindel file');

my $output_dir = Genome::Sys->create_temp_directory();

my $cmd = $class->create(
    input_dindel_file => $input_file,
    output_directory => $output_dir,
    num_windows_per_file => -1,
);
ok($cmd->execute(), "Successfully ran command");

my $expected_output_dir = File::Spec->join($test_dir, 'expected-output');
test_txt_files_are_identical();

done_testing();

sub test_txt_files_are_identical{
    my @expected_files = sort(get_text_files($expected_output_dir));

    for my $file (@expected_files) {
        my $found_file = $file;
        $found_file =~ s/$expected_output_dir/$output_dir/;
        compare_ok($file, $found_file, "$found_file identical to $file");
    }
}

sub get_text_files {
    my $base = shift;
    my @files = glob(File::Spec->join($base, '*.txt'));
    push @files, glob(File::Spec->join($base, '*' ,'*.txt'));
    push @files, glob(File::Spec->join($base, '*', '*', '*.txt'));
    return @files;
}
