#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
require File::Compare;

use_ok( 'Genome::Model::Tools::Consed::RemoveTags' ) or die;

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Consed-RemoveTags';
ok( -d $data_dir, 'test suite dir exists') or die;

my $ace_in = '454Contigs.ace';
my $ace_out = $ace_in.'.tags_removed';
ok( -s $data_dir.'/'.$ace_in, 'Test ace file exists' ) or die;

my $test_dir = Genome::Sys->create_temp_directory();
ok( -d $test_dir, 'made temp test dir' );

my @tags = qw/ tear join /;

my $create = Genome::Model::Tools::Consed::RemoveTags->create(
    ace_in => $data_dir.'/'.$ace_in,
    ace_out => $test_dir.'/'.$ace_out,
    tags => \@tags,
);
ok( $create, 'created tool' );
ok( $create->execute, 'successfully executed tool' );

ok( File::Compare::compare($data_dir.'/'.$ace_out,$test_dir.'/'.$ace_out)== 0, 'Output files match' );

#<STDIN>;

done_testing();

exit;
