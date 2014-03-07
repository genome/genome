#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Test::Exception;
use Genome::Utility::Test 'compare_ok';
use File::Spec;
use Sub::Install qw();
use File::Basename qw(dirname);

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

class ParkTest {
    is => 'Genome::Model::Tools::Park::Base',
    has => [
        _template_path => {
            is => 'Path',
        },
        bar => {
            is => 'Number',
        },
        bof => {
            is => 'Number',
        },
        blah => {
            is => 'Number',
        },
        blerg => {
            is => 'Number',
            is_many => 1,
        },
    ],
};

my $test_dir = __FILE__ . '.d';
my $test_template = File::Spec->join($test_dir, 'test_template.tsv');
my $bad_template = File::Spec->join($test_dir, 'bad_template.tsv');

test_basic_usage();
test_missing_accessor();
test_undef_attribute();

test_run_generic_process();
test_run_rex_process();

test_rex_script_path();

done_testing();


sub test_basic_usage {
    my $test_class = ParkTest->create(
        _template_path => $test_template,
        bar => 1,
        bof => 2,
        blah => 3,
        blerg => [4,5],
    );
    my $inputs_filename = $test_class->_generate_inputs_file();

    my $expected_inputs = File::Spec->join($test_dir, 'expected_inputs.tsv');
    compare_ok($inputs_filename, $expected_inputs, 'generated the expected inputs file');
}

sub test_missing_accessor {
    my $bad_class = ParkTest->create(
        _template_path => $bad_template,
    );
    dies_ok {$bad_class->_generate_inputs_file();} "dies as expected if accessor is missing";
    note $@;
}

sub test_undef_attribute {
    my $undef_class = ParkTest->create(
        _template_path => $test_template,
        bar => 1,
        bof => 2,
        blerg => [4,5],
    );
    dies_ok {$undef_class->_generate_inputs_file();} "dies as expected if attribute is undef";
    note $@;
}

sub test_rex_script_path {
    my $test_class = ParkTest->create();
    my $result = $test_class->_rex_script_path('process start');
    like($result, qr(rex_process_start\.sh$), "got expected script ($result) from command ('process_start')");
    dies_ok {$test_class->_rex_script_path('missing command');} "fails to find missing command";
    note $@;
}

sub test_run_generic_process {
    my $test_class = ParkTest->create();
    my $process_uri = '/this/is/a/test/';
    my $dir = dirname(__FILE__);
    my $cmd = [File::Spec->join($dir, 'test.sh')];
    is($test_class->_run_generic_process($cmd), $process_uri,
        "Got process uri ($process_uri) from run_generic_process");
}

sub test_run_rex_process {
    Sub::Install::reinstall_sub({
        into => 'Genome::Model::Tools::Park::Base',
        as => '_run_generic_process',
        code => sub { my $self = shift; return shift; },
    });
    my $test_class = ParkTest->create();
    my $source_path = 'Test::Path';
    my $inputs_filename = '/test/filename.tsv';
    my $cmd = $test_class->_run_rex_process($source_path, $inputs_filename);
    is(scalar(@$cmd), 5, "found expected number of elements in \$cmd");
    shift @$cmd;
    is_deeply($cmd, ['--source-path', $source_path, '--inputs-file', $inputs_filename], 'found expected \$cmd');
}
