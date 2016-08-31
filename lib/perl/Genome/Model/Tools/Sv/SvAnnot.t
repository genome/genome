#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 7;
use File::Compare;
use File::Temp;

use_ok('Genome::Model::Tools::Sv::SvAnnot');

my $test_input_dir  = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-Sv-SvAnnot/v2/';
my $tmp_dir = Genome::Sys->create_temp_directory();

for my $type (qw(36_standard 37_merged)) {
    my $sv_file    = $test_input_dir . 'sv.file.' . $type;
    my $expect_out = $test_input_dir . 'sv.annot.'. $type;
    my $out_file   = $tmp_dir . '/sv.annot.' . $type;

    my ($build, $format) = $type =~ /^(\d+)_(\S+)$/;
    my %params = (
        sv_file     => $sv_file,
        sv_format   => $format,
        output_file => $out_file,
        annot_build => $build,
    );

    my $annot_valid = Genome::Model::Tools::Sv::SvAnnot->create(%params);
   
    ok($annot_valid, "created SvAnnot object ok for $type");
    ok($annot_valid->execute(), "executed SvAnnot object OK for $type");
    is(compare($out_file, $expect_out), 0, "output matched expected result for $type");
}
