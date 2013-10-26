#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 9;
use File::Compare;

use above 'Genome';

BEGIN {
    use_ok('Genome::Model::Tools::Fastq::Trimq2::Simple');
};

my $base_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fastq/Trimq2/Simple';

my $tmp_dir = File::Temp::tempdir(
    'Fastq-Trimq2-Simple-XXXXX', 
    TMPDIR => 1,
    CLEANUP => 1
);
my $fastq_file = "$base_dir/test_simple.fastq";

#Test hard trim
my $trimq2 = Genome::Model::Tools::Fastq::Trimq2::Simple->create(
    fastq_file  => $fastq_file,
    output_dir  => $tmp_dir,
    trim_style  => 'hard',
    #trim_string => '#',  default is #
);
isa_ok($trimq2,'Genome::Model::Tools::Fastq::Trimq2::Simple');

ok($trimq2->execute,'execute command '. $trimq2->command_name);

is($trimq2->out_file, "$tmp_dir/test_simple.trimq2.fastq", 'output name is ok');

compare_output('hard');

#Test Bwa trim

$trimq2 = Genome::Model::Tools::Fastq::Trimq2::Simple->create(
    fastq_file  => $fastq_file,
    output_dir  => $tmp_dir,
    trim_style  => 'smart1',
    #trim_string => '#',  default is #
);
ok($trimq2->execute,'execute command '. $trimq2->command_name);

compare_output('smart1');

sub compare_output {
    my $style = shift;

    for my $file qw(test_simple.trimq2.fastq trimq2.report) {
        my $output_file = $tmp_dir."/$file";
        my $expect_file = $base_dir."/$file.$style";
        ok(compare($output_file, $expect_file) == 0, "Output $file with $style style is created as expected");
        unlink $output_file if $style eq 'hard'; #make file name available for next round
    }
}
