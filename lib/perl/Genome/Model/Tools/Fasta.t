#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Cwd 'getcwd';
use Data::Dumper 'Dumper';
use Test::More tests => 39;
use File::Temp 'tempdir';
use File::Compare 'compare';
use File::Copy 'copy';
use Bio::SeqIO;

use_ok ('Genome::Model::Tools::Fasta')
    or die;

#< SETUP (8) >#
my $CWD = getcwd();
ok(-d $CWD, "We are in dir: $CWD");
# Original files
my $ORIG_DIR = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-Fasta';
ok(-d $ORIG_DIR, "Test dir ($ORIG_DIR) exists");
my $ORIG_FASTA = $ORIG_DIR.'/file.fasta';
ok(-f $ORIG_FASTA, "Original fasta ($ORIG_FASTA) exists");
my $ORIG_QUAL = $ORIG_FASTA.'.qual';
ok(-f $ORIG_QUAL, "Original qual ($ORIG_QUAL) exists");
# Temp files 
# not using a tmp prefix for convenience
my $DIR = tempdir(CLEANUP => 1);
ok(-d $DIR, "Temp dir ($DIR) exists");
my $FASTA = $DIR.'/file.fasta';
File::Copy::copy($ORIG_FASTA, $FASTA);
ok(-f $FASTA, "Temp fasta ($FASTA) exists");
my $QUAL = $FASTA.'.qual';
File::Copy::copy($ORIG_QUAL, $QUAL);
ok(-f $QUAL, "Temp qual ($QUAL) exists");

my $fasta_command = Genome::Model::Tools::Fasta->create(
    fasta_file => $FASTA,
);
ok($fasta_command, "Created - let's go!")
    or die;

#< QUAL FILE (1) >#
ok($fasta_command->have_qual_file, 'Got the quality file');

#< HELPS (2) >#
ok($fasta_command->help_brief, "Got helped, in brief");
ok($fasta_command->help_detail, "Got helped, in detail");

#< CHDIR METHODS (2) >#
$fasta_command->chdir_cwd;
is(getcwd(), $CWD, "Changed back to dir wher etest was launched");
$fasta_command->chdir_fasta_directory;
is(getcwd(), $DIR, "Changed to the tmp dir for testing");

#< FILE METHODS (6) >#
my %file_tests = (
    fasta_directory => $DIR.'/',
    fasta_basename => 'file',
    fasta_suffix => '.fasta',
    fasta_base => 'file.fasta',
    qual_base => 'file.fasta.qual',
    qual_file => $FASTA.'.qual',
);
for my $test ( keys %file_tests ) {
    is($fasta_command->$test, $file_tests{$test}, "Matched: $test");
}

#< SUFFIX (4) >#
my $suffix = 'sanitized';
my %suffix_tests = (
    fasta_base_with_new_suffix => "file.$suffix.fasta",
    fasta_file_with_new_suffix => $DIR."/file.$suffix.fasta",
    qual_base_with_new_suffix => "file.$suffix.fasta.qual",
    qual_file_with_new_suffix => $DIR."/file.$suffix.fasta.qual",
);
for my $test ( keys %suffix_tests ) {
    is($fasta_command->$test($suffix), $suffix_tests{$test}, "Matched: $test");
}

#< BACK UP FILES (8) >#
# Names
is($fasta_command->default_back_up_suffix, 'bak', "Back up suffix");
my $back_suffix = 'back';
my %back_up_name_tests = (
    fasta_back_up_file => [qw/ file.bak.fasta file.back.fasta /],
    qual_back_up_file =>  [qw/ file.bak.fasta.qual file.back.fasta.qual /],
);
for my $test ( keys %back_up_name_tests ) {
    is($fasta_command->$test, $DIR.'/'.$back_up_name_tests{$test}->[0], "Matched: $test");
    is($fasta_command->$test($back_suffix), $DIR.'/'.$back_up_name_tests{$test}->[1], "Matched: $test");
}
# Do a back up
my @bak_files = $fasta_command->back_up_fasta_and_qual_files;
ok(@bak_files, 'Backed up fasta and qual files');
is(compare($FASTA, $bak_files[0]), 0, 'Backed up fasta matches');
is(compare($QUAL, $bak_files[1]), 0, 'Backed up qual matches');

#< BIOSEQ (4) >#
my %bioseq_tests = (
    get_fasta_reader => $FASTA,
    get_qual_reader => $QUAL,
    get_fasta_writer => $DIR.'/new.fasta',
    get_qual_writer => $DIR.'/new.fasta.qual',
);
for my $test ( keys %bioseq_tests ) {
    isa_ok($fasta_command->$test($bioseq_tests{$test}), 'Bio::SeqIO', "Got Bio::SeqIO for $test");
}

#< CREATE FAILS >#
ok( ! Genome::Model::Tools::Fasta->create(), "Create w/o fasta file - failed as expected");
ok( ! Genome::Model::Tools::Fasta->create(fasta_file => $ORIG_DIR.'/file.txt'), "Create w/ non existing fasta - failed as expected");
ok( ! Genome::Model::Tools::Fasta->create(fasta_file => $ORIG_DIR.'/file.fake'), "Create w/ invalid fasta ext - failed as expected");
chdir $CWD;
