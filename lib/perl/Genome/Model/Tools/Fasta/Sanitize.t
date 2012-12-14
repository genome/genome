#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Bio::SeqIO;
require Cwd;
use Data::Dumper 'Dumper';
use File::Compare 'compare';
use File::Temp 'tempdir';
use File::Basename;
use Test::More tests => 30;

#< README >#
# All of the expected fasta files were made by hand

use_ok ('Genome::Model::Tools::Fasta::Sanitize')
    or die;

##< SETUP: 3 TESTS >#
my $DIR = Cwd::abs_path($ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fasta-Sanitize');
ok(-d $DIR, "Test dir ($DIR) exists");
my $DIRTY_FASTA = $DIR.'/dirty.fasta';
ok(-f $DIRTY_FASTA, "Dirty fasta ($DIRTY_FASTA) exists");
#my $TMPDIR = tempdir(DIR => $DIR, CLEANUP => 1);
my $TMPDIR = Genome::Sys->create_temp_directory;
ok(-d $TMPDIR, "Temp dir ($TMPDIR) exists");

my $san = Genome::Model::Tools::Fasta::Sanitize->create(
    fasta_file => $DIRTY_FASTA,
    capitalize => 1, # gotta include sumpin'
);
ok($san, "Created sanitize for testing default sanitized fasta name")
    or die;

#< MISC (4) >#
ok($san->help_brief, "Got helped, in brief");
ok($san->help_detail, "Got helped, in detail");
ok($san->sanitation_methods, 'Got sanitation methods');
is($san->default_sanitized_fasta_file, $DIR.'/dirty.sanitized.fasta', 'Default sanitized fasta file name is correct');

#< SANITIZERS, SINLGY >#
my %sanitizers_and_params = (
    capitalize => 1,
    min_threshold => 1000,
    max_threshold => 1000,
    replace_xs_with_ns => 1,
    replace => 'a{22,}:'.('NX' x 11),
    rm_descs => 1,
);

for my $sanitizer ( keys %sanitizers_and_params ) {
    my $sanitized_fasta = sprintf('%s/sanitized.%s.fasta', $TMPDIR, $sanitizer);
    my $expected_fasta = sprintf('%s/expected.%s.fasta', $DIR, $sanitizer);
    $san = Genome::Model::Tools::Fasta::Sanitize->create(
        fasta_file => $DIRTY_FASTA,
        sanitized_fasta_file => $sanitized_fasta,
        $sanitizer => $sanitizers_and_params{$sanitizer},
    );
    ok($san, "Created sanitize for $sanitizer")
        or die;
    ok($san->execute, "Executed");
    is(compare($sanitized_fasta, $expected_fasta), 0, 'Expected and generated fasta matches');
}

#< CREATE FAILS >#
ok( ! Genome::Model::Tools::Fasta::Sanitize->create(fasta_file => $DIRTY_FASTA), "Create w/ no sanitation to do - failed as expected");
ok( 
    ! Genome::Model::Tools::Fasta::Sanitize->create(fasta_file => $DIRTY_FASTA, replace => 'no_replace:'),
    "Create w/ invalid replace string - failed as expected"
);
ok( 
    ! Genome::Model::Tools::Fasta::Sanitize->create(fasta_file => $DIRTY_FASTA, replace => ':no_find'),
    "Create w/ invalid replace string - failed as expected"
);



#$HeadURL$
#$Id$
