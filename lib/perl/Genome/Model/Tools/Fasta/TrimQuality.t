#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 11;
use File::Compare;
use File::Temp;
use File::Copy;

BEGIN {
        use_ok ('Genome::Model::Tools::Fasta::TrimQuality');
}

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fasta/TrimQuality';
ok(-d $dir, "Test dir ($dir) exists");
my $fasta = $dir .'/test.fasta';
ok(-f $fasta, "Fasta ($fasta) exists");
my $qual = $fasta .'.qual';
ok(-f $qual, "Qual ($qual) exists");

my $tmp_dir  = File::Temp::tempdir(
    "Fasta_TrimQuality_XXXXXX", 
    DIR     => "$ENV{GENOME_TEST_TEMP}",
    CLEANUP => 1,
);

my @test_tmp_dirs = map{$tmp_dir.'/test'.$_}qw(1 4 5);
map{mkdir $_}@test_tmp_dirs;

my @test_fastas;
my @test_quals;

for my $test_tmp_dir (@test_tmp_dirs) {
    copy $dir.'/test.fasta', $test_tmp_dir.'/test.fasta';
    push @test_fastas, $test_tmp_dir.'/test.fasta';
    copy $dir.'/test.fasta.qual', $test_tmp_dir.'/test.fasta.qual';
    push @test_quals, $test_tmp_dir.'/test.fasta.qual';
}   

# Should work
my $trim1 = Genome::Model::Tools::Fasta::TrimQuality->create(
    fasta_file => $test_fastas[0],
    min_trim_quality => 10,
    min_trim_length  => 100,
);
ok($trim1->execute, "trim1 finished ok");

# No fasta
ok(
    ! Genome::Model::Tools::Fasta::TrimQuality->create(
        fasta_file => $dir.'/no.fasta',
        min_trim_quality => 12,
        min_trim_length  => 80,
    ),
    "This supposed to fail because of no fasta",
);

# No qual for fasta 
ok(
    ! Genome::Model::Tools::Fasta::TrimQuality->create(
        fasta_file => $dir.'/test.no_qual.fasta',
        min_trim_quality => 12,
        min_trim_length  => 80,
    ),
    "This supposed to fail because of no quality file",
);

my $trim4 = Genome::Model::Tools::Fasta::TrimQuality->create(
    fasta_file => $test_fastas[1],
);

ok($trim4->execute, "trim4 running ok");
is(compare("$dir/test.fasta.ori.clip", $test_fastas[1]),0, "using default, fasta is ok");
cmp_ok(compare("$dir/test.fasta.qual.ori.clip", $test_quals[1]),'==', 0, "using default, qual is ok");


my $trim5 = Genome::Model::Tools::Fasta::TrimQuality->create(
    fasta_file       => $test_fastas[2],
    min_trim_length => 'hello',
);

#eval(my $rv = $trim5->execute);
#print $rv;
ok(!$trim5->execute, "This supposed to fail because trim5 uses non-integer as min_trim_length");

exit;


