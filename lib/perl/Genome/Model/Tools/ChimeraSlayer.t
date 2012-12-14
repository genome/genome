#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;
require File::Compare;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

use_ok( 'Genome::Model::Tools::ChimeraSlayer' );

my $version = 1;
my $test_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ChimeraSlayer/v'.$version;
ok( -d $test_data_dir, 'Test data dir' );

#check test files
my $test_input = 'chims1.NAST';
ok( -s $test_data_dir.'/'.$test_input, 'Test input file' );

my @test_outputs = qw/
chims1.NAST.CPS
chims1.NAST.CPS.CPC
chims1.NAST.CPS.CPC.align
chims1.NAST.CPS.CPC.wTaxons
chims1.NAST.CPS_RENAST
chims1.NAST.CPS_RENAST.cidx
chims1.NAST.cidx
/;
for my $file ( @test_outputs ) {
    ok ( $test_data_dir.'/'.$file, "Test output file: $file" );
}

#copy test file
my $temp_test_dir = Genome::Sys->create_temp_directory();
ok( -d $temp_test_dir, 'Temp test dir' );
ok( File::Copy::copy( $test_data_dir.'/'.$test_input, $temp_test_dir ), 'Copied test input file' );

#create/execute tool
my $tool = Genome::Model::Tools::ChimeraSlayer->create(
    query_NAST => $temp_test_dir.'/'.$test_input,
    exec_dir => $temp_test_dir,
    printCSalignments => 1,
    printFinalAlignments => 1,
);
ok( $tool, 'Created tool' );
ok( $tool->execute, 'Executed tool' );

#check outputs
#files that should consistantly match
for my $file ( qw/ chims1.NAST.CPS chims1.NAST.CPS_RENAST / ) {
    ok( -s $temp_test_dir.'/'.$file, "Created file: $file" );
    ok( File::Compare::compare( $temp_test_dir.'/'.$file, $test_data_dir.'/'.$file ) == 0, "Test $file files match" );
}
#files that differ sometimes
 #chims1.NAST.CPS.CPC.align
 #chims1.NAST.CPS.CPC
 #chims1.NAST.CPS.CPC.wTaxons
 #chims1.NAST.CPS_RENAST
#binary files
 #chims1.NAST.CPS_RENAST.cidx
 #chims1.NAST.cidx

#<STDIN>;

done_testing();
