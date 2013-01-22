#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Temp;
use Test::More;# tests => 4;
use Genome::Utility::Test qw(compare_ok);

use_ok('Genome::Model::Tools::Velvet::ToAce');

#test suite dir
my $test_version = '-v3';
my $root_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Velvet/ToAce'.$test_version;

#test dir
my $run_dir = Genome::Sys->create_temp_directory();
ok( -d $run_dir, "Create test dir" );

#link input files
for my $file (qw/  Sequences velvet_asm.afg / ) {
    ok( -s $root_dir."/$file", "input $file exists" );
    symlink( $root_dir."/$file", $run_dir."/$file" );
    ok( -s $root_dir."/$file", "linked input file" );
}

#create/execute tool
my $ta = Genome::Model::Tools::Velvet::ToAce->create(
    assembly_directory => $run_dir,
    time        => 'Wed Jul 29 10:59:26 2009',
    #sqlite_yes => 1,
);

ok($ta, 'to-ace creates ok');
ok($ta->execute, 'velvet to-ace runs ok');

#check ace file
my $out_ace = $run_dir.'/edit_dir/velvet_asm.ace';
ok( -s $out_ace, "Created ace file" );

my $ori_ace = $root_dir.'/velvet_asm.ace';
ok( -s $ori_ace, "Test ace file exists" );

compare_ok($ori_ace, $out_ace, 'Ace file converted from velvet output is OK',
    filters => [qr(^comment VelvetToAce.*), qr(^Run by .*)]);

#<STDIN>;

done_testing();
