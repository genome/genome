#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 6;
use File::Temp;

BEGIN {
        use_ok('Genome::Model::Tools::BioDbFasta::Convert');
}

my $infile_nonexist = "/blah/not/exists";
#my $outfile = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-BioDbFasta/Convert/testoutput1";

my $tmp_dir = File::Temp::tempdir(
    "BioDbFasta_convert_XXXXXX", 
    TMPDIR => 1,
    CLEANUP => 1,
);

my $outfile = $tmp_dir.'/testoutput1';

my $infile_test = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-BioDbFasta/Convert/testinput1";

my $c_f_noexist = Genome::Model::Tools::BioDbFasta::Convert->create('infile' => $infile_nonexist  ,
                                                          'outfile' =>  $outfile );

isa_ok($c_f_noexist,'Genome::Model::Tools::BioDbFasta::Convert');

is($c_f_noexist->execute(),0,'infile non-existent');

my $c =  Genome::Model::Tools::BioDbFasta::Convert->create('infile' => $infile_test ,
                                                          'outfile' =>  $outfile );

ok($c->execute(),'running on test file 1');

my $c_test_gqs = Genome::Model::Tools::BioDbFasta::Convert->create('infile' => $infile_test ,
                                                          'outfile' =>  $outfile );

my $string = $c_test_gqs->getQualString("55 33 55 55 10 29 34");
is($string,'XBXX+=C','getQualString');


my $c_out_no_x = Genome::Model::Tools::BioDbFasta::Convert->create('infile' => $infile_test ,
                                                          'outfile' =>  '/' );

my $ret;
eval { $ret = $c_out_no_x->execute(); };
is($ret, undef, 'trying to open file that can\'t be opened');
