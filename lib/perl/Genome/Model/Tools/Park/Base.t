#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Test::Exception;
use Genome::Utility::Test 'compare_ok';
use File::Spec;

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
    compare_ok($inputs_filename, $expected_inputs, 'Generated the expected inputs file');
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
