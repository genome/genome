#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 12;
use Genome::Utility::Test;
use File::Spec;
use File::Compare;

my $class ='Genome::Model::Tools::CopyCat::Somatic'; 
use_ok($class);

my $test_version = 1;
my $test_dir = File::Spec->join(Genome::Utility::Test->data_dir($class), "v$test_version");
ok(-e $test_dir, "Test directory $test_dir exists");

my $tumor_window_file = File::Spec->join($test_dir, 'tumor.wind');
ok(-s $tumor_window_file, "tumor wind file exists");

my $normal_window_file = File::Spec->join($test_dir, 'normal.wind');
ok(-s $normal_window_file, "normal wind file exists");

my $tumor_samtools_file = File::Spec->join($test_dir, 'tumor_samtools', 'snvs.hq');
ok(-s $tumor_samtools_file, 'tumor samtools file exists');

my $normal_samtools_file = File::Spec->join($test_dir, 'normal_samtools', 'snvs.hq');
ok(-s $normal_samtools_file, 'normal samtools file exists');

my $output_directory = Genome::Sys->create_temp_directory();
my $annotation_directory = '/gscmnt/gc6122/info/medseq/annotations/copyCat'; #TODO: don't do this, this is really bad

my $cmd = Genome::Model::Tools::CopyCat::Somatic->create(
    normal_window_file => $normal_window_file,
    tumor_window_file => $tumor_window_file,
    output_directory => $output_directory,
    annotation_directory => $annotation_directory,
    per_library => 1,
    per_read_length => 1,
    genome_build => 'hg19.chr14only',
    normal_samtools_file => $normal_samtools_file,
    tumor_samtools_file => $tumor_samtools_file,
);
ok($cmd, "Created command successfully");
ok($cmd->execute, "Executed the command successfully");

my $expected_output_dir = File::Spec->join($test_dir, 'expected');
my @diffable_files = qw| libraries.tumor.txt
                         libraries.normal.txt |;

for my $file (@diffable_files){
    my $expected = File::Spec->join($expected_output_dir, $file);
    my $actual = File::Spec->join($output_directory, $file);
    is(compare($actual, $expected),0,"Actual file is the same as the expected file: $file");
}

my @non_diffable_files = qw| alts.paired.dat
                             segs.paired.dat |;

for my $file (@non_diffable_files){
    my $expected = File::Spec->join($expected_output_dir, $file);
    my $actual = File::Spec->join($output_directory, $file);
    my ($actual_wc) = split(" ", `wc -l $actual`);    
    my ($expected_wc) = split(" ", `wc -l $expected`);    
    ok(abs ($expected_wc - $actual_wc) <= 1, "$file line length is withing tolerance");
}

1;
