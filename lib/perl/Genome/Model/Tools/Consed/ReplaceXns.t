#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::Consed::ReplaceXns' ) or die;

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Consed-ReplaceXns';
ok( -d $data_dir, 'test suite dir' );

my $ace_in = '454Contigs.ace.1';
my $ace_out = '454Contigs.ace.1.xns_replaced';

ok( -s $data_dir.'/'.$ace_in, 'test suite ace in' );
ok( -s $data_dir.'/'.$ace_out, 'test suite ace out' );

my $test_dir = Genome::Sys->create_temp_directory();
ok( -d $test_dir, 'temp test dir' );

my $create = Genome::Model::Tools::Consed::ReplaceXns->create(
    ace_in  => $data_dir.'/'.$ace_in,
    ace_out => $test_dir.'/'.$ace_out,
);
ok( $create, 'created tool' );
ok( $create->execute, 'executed tool' );

ok( -s $test_dir.'/'.$ace_out, 'created ace out' );

#output file will vary because a base can be selected randomly
#when all reads at consensus call Ns or when there more than one
#top base calls .. there are 11 instances in test ace where a base
#is selected randomly

my $diff_command = 'sdiff -s '.$data_dir.'/'.$ace_out.' '.$test_dir.'/'.$ace_out;
my @ace_diffs = `$diff_command`;
for my $diff ( @ace_diffs ) {
    chomp $diff;
    #print $diff."\n";
    ok ( $diff =~ /^[acgtxn*\s\|]+$/i, 'diff is due to base differences caused by random base selection' );
}

#<STDIN>;

done_testing();
