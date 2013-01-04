#!/usr/bin/env genome-perl

use strict;
use warnings FATAL => 'all';

use above 'Genome';
use Test::More;

sub get_cmd_string_for_hash {
    my $params = shift;
    my $cmd = Genome::Model::Tools::Picard::FilterSamReads->create(%{$params});
    return $cmd->_generate_cmd_string("");
}

my $params1 = {
    filter      => 'includeAligned',
    input_file  => 'input_test',
    output_file => 'output_test',
};
my $cmd1_expected = " net.sf.picard.sam.FilterSamReads " .
    "O=output_test I=input_test " .
    "FILTER=includeAligned " .
    "WRITE_READS_FILES=true";
my $cmd1_actual = get_cmd_string_for_hash( $params1 );
is( $cmd1_expected,  $cmd1_actual, "Expected command matches output" );

my $params2 = {
    filter      => 'includeAligned',
    input_file  => 'input_test',
    output_file => 'output_test',
    sort_order  => 'unsorted',
    use_version => '1.77',
};
my $cmd2_expected = " net.sf.picard.sam.FilterSamReads " .
    "O=output_test I=input_test " .
    "FILTER=includeAligned " .
    "WRITE_READS_FILES=true " .
    "SORT_ORDER=unsorted";
my $cmd2_actual = get_cmd_string_for_hash( $params2 );
is( $cmd2_expected, $cmd2_actual, "Expected command matches output" );

eval {
    my $filter_cmd = Genome::Model::Tools::Picard::FilterSamReads->create(
        filter      => 'includeAligned',
        input_file  => 'input_test',
        output_file => 'output_test',
        sort_order  => 'unsorted',
    );

    $filter_cmd->execute();
    die;
};

if ( $@ ) {
    my $error_str = $@;
    ok( $error_str =~ m/may not be available/,
      "Crashed with wrong Picard version");
}

done_testing();

1;
