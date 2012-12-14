#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Compare;
use Test::More;

use_ok('Genome::Model::Tools::Sx::Trim::Quake') or die;

my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/';
my $input = $test_dir.'/reader_writer.collated.fastq';
my $tmp_dir = Genome::Sys->base_temp_directory;
my $output = $tmp_dir.'/output.fastq';

my %quake_params = map { $_ => 1 } Genome::Model::Tools::Sx::Trim::Quake->quake_param_names;
$quake_params{input} = [$input.':cnt=2'];
$quake_params{output} = [$output];
my $quake;
no warnings;
*Genome::Sys::shellcmd = sub{ 
    # DOES NOT RUN QUAKE!
    my ($self, %params) = @_;
    is(
        $params{cmd},
        'quake.py -q 33 -r '.$quake->_tmpdir.'/quake.fastq --hash_size 1 --headers --int -k 1 -l 1 --log --no_count --no_cut --no_jelly -p 1 --ratio 1 -t 1 -u',
        'quake command matches',
    );
    Genome::Sys->copy_file($input, $quake->_tmpdir.'/quake.cor.fastq');
    return 1; 
};
use warnings;

$quake = Genome::Model::Tools::Sx::Trim::Quake->create(
    %quake_params,
);
ok($quake, 'create');
$quake->dump_status_messages(1);
ok($quake->execute, 'execute');
is(File::Compare::compare($input, $quake->_tmpdir.'/quake.fastq'), 0, 'input written correctly');
is(File::Compare::compare($output, $input), 0, 'output file matches');

#print "$tmp_dir\n";<STDIN>;
done_testing();
