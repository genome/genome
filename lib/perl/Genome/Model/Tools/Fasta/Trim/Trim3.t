#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use File::Compare 'compare';
use File::Temp 'tempdir';
use File::Copy 'copy';
use Test::More;

use_ok ('Genome::Model::Tools::Fasta::Trim::Trim3') or die;

my $dir = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-Fasta/TrimQuality';
ok(-d $dir, "Test dir ($dir) exists");
my $example_fasta = $dir .'/example.fasta';
ok(-f $example_fasta, "Example fasta ($example_fasta) exists");
my $example_qual = $example_fasta .'.qual';
ok(-f $example_qual, "Example qual ($example_qual) exists");
my $result_fasta = $dir .'/result.fasta';
ok(-f $result_fasta, "result fasta ($result_fasta) exists");
my $result_qual = $result_fasta .'.qual';
ok(-f $result_qual, "result qual ($result_qual) exists");

my $tmpdir = tempdir(CLEANUP => 1);
my $fasta = $tmpdir.'/example.fasta';
my $qual = $fasta.'.qual';

# Fails - no fasta
ok(
    !Genome::Model::Tools::Fasta::Trim::Trim3->create(
        fasta_file => $fasta,
        min_trim_quality => 12,
        min_trim_length => 80,
    ),
    "Failed as expected - fasta doesn't exist",
);

# Copy fasta
Genome::Sys->copy_file($example_fasta, $fasta);
ok(-s $fasta, 'Copied fasta');

# Fails - no qual
ok(
    !Genome::Model::Tools::Fasta::Trim::Trim3->create(
        fasta_file => $fasta,
        min_trim_quality => 12,
        min_trim_length => 80,
    ),
    "Failed as expected - no qual file",
);

# Copy qual
Genome::Sys->copy_file($example_qual, $qual);
ok(-s $qual, 'Copied qual');

# Fail invalid params
ok(
    !Genome::Model::Tools::Fasta::Trim::Trim3->create(
        fasta_file => $fasta,
        min_trim_length => 'hello',
    ),
    'Failed as expected invalid min_trim_length',
);
ok(
    !Genome::Model::Tools::Fasta::Trim::Trim3->create(
        fasta_file => $fasta,
        min_trim_length => 'hello',
    ),
    'Failed as expected invalid min_trim_length',
);

# Valid (default params)
my $trim = Genome::Model::Tools::Fasta::Trim::Trim3->create(
    fasta_file => $fasta,
);
ok($trim, "Create");
ok($trim->execute, "Execute");
is(compare($fasta, $result_fasta), 0, 'Fasta matches');
is(compare($qual, $result_qual), 0, 'Qual matches');

#print "$tmpdir\n"; <STDIN>;
done_testing();
