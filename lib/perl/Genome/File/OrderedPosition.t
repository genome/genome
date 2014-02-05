#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw/compare_ok/;

use strict;
use warnings;

my $pkg = 'Genome::File::OrderedPosition';

use_ok($pkg);

my $TEST_VERSION = 2;

# my $data_dir_t = Genome::Utility::Test->data_dir_ok($pkg, $TEST_VERSION);
my $data_dir = File::Spec->join($ENV{GENOME_TEST_INPUTS}, "Genome-File-OrderedPosition", $TEST_VERSION);

subtest "sorted file without header" => sub {
    my $file = $pkg->new(File::Spec->join($data_dir, "varscan.snp.Somatic.strandfilter"), 10_000 );
    ok(!defined($file->{header}), "We didn't find a header");
    is($file->{line_number}, 1, "We're on the first line");

    my $line = $file->getline_for_position(1, 177);
    ok(defined($line), "We found a line with position (1, 177)");
    is($file->{line_number}, 1, "We found the position (1, 177) on the first line");
    is(
        $line,
        "1\t177\tA\tC\t228\t24\t9.52%\tM\t240\t37\t13.36%\tM\tGermline\t6.917202812301296E-20\t0.1067345298531505\t144\t96\t25\t12\t0.40\t0.49\t0.60\t0.68\t13.80\t14.22\t0.42\t0.52\t0.52\t44.28\t47.57\n",
        "First line as expected"
    );

    $line = $file->getline_for_position(1, 177);
    ok(defined($line), "We found the same line again: position (1, 177)");

    $line = $file->getline();
    ok(defined($line), "We got a line");
    is($file->{line_number}, 2, "We're on the second line");
    is(
        $line,
        "1\t180\tT\tC\t350\t36\t9.33%\tY\t283\t43\t13.19%\tY\tSomatic\t1.6733783183436353E-25\t0.06501933029082527\t178\t105\t18\t25\t0.35\t0.48\t0.63\t0.42\t9.24\t19.28\t10.04\t0.58\t0.51\t44.73\t47.23\n",
        "Second line as expected"
    );

    $line = $file->getline_for_position(1, 327);
    ok(defined($line), "We found a line with position (1, 327)");
    is($file->{line_number}, 5, "We're on the fifth line");
    is(
        $line,
        "1\t327\tT\tC\t30\t11\t26.83%\tY\t83\t26\t23.85%\tY\tGermline\t5.763307701930354E-13\t0.7255614827940279\t37\t46\t9\t17\t0.48\t0.58\t0.45\t0.35\t22.01\t37.58\t15.57\t0.60\t0.50\t41.59\t45.46\n",
        "Fifth line as expected"
    );

    $line = $file->getline_for_position(1, 177);
    ok(!defined($line), "Line not found (previous position): position (1,177)");

    $line = $file->getline_for_position(1, 444);
    ok(!defined($line), "Line not found (position non-existent): position (1, 444)");

    $line = $file->getline_for_position(1, 998);
    ok(!defined($line), "Line not found (position non-existent after EOF): position (1, 998)");
};

subtest "sorted file with header" => sub {
    my $file_with_header = $pkg->new(File::Spec->join($data_dir, "varscan.snp"), 10_000 );
    ok(defined($file_with_header->{header}), "We found a header");

    my $line = $file_with_header->getline_for_position(1, 177);
    ok(defined($line), "We found a line with position (1, 177)");
    is($file_with_header->{line_number}, 176, "We're on the 176th line");
    is(
        $line,
        "1\t177\tA\tC\t228\t24\t9.52%\tM\t240\t37\t13.36%\tM\tGermline\t6.917202812301296E-20\t0.1067345298531505\t144\t96\t25\t12\n",
        "176th line as expected"
    );

    $line = $file_with_header->getline_for_position("X", 2);
    ok(defined($line), "We found a line with position ('X', 2)");
    is($file_with_header->{line_number}, 707, "We're on the 707th line");
};

subtest "unsorted file with header" => sub {
    my $unsorted_file = $pkg->new(File::Spec->join($data_dir, "varscan.snp.unsorted"), 10_000 );
    compare_ok($unsorted_file->{sorted_filename}, File::Spec->join($data_dir, "varscan.snp.without_header"), "Sorting works correctly");
};

done_testing();
