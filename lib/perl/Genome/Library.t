#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above "Genome";

use Test::More tests => 3;

subtest 'use and roles' => sub{
    plan tests => 2;

    use_ok('Genome::Library') or die;
    ok(Genome::Library->does('Genome::Role::Notable'), 'Library does Genome::Role::Notable');

};

my $library;
subtest 'create' => sub{
    plan tests => 4;

    my $sample = Genome::Sample->create(name => '__TEST_SAMPLE__');
    ok($sample, 'create sample');

    $library = Genome::Library->create(
        sample => $sample,
        name => $sample->name . "-extlibs",
        original_insert_size => '1kb',
        library_insert_size => '300-500',
        protocol => 'karate chop',
    );
    isa_ok($library, 'Genome::Library');
    is($library->name, $sample->name . "-extlibs", "name is what is expected");

    my $commit = eval{ UR::Context->commit; };
    ok($commit, 'commit');

};

subtest 'is_rna' => sub{
    plan tests => 4;

    $library->sample->extraction_type('dna');
    is($library->is_rna, 0, 'is NOT rna when extraction_type is dna and transcript_strand is UNDEF');

    $library->transcript_strand('unstranded');
    is($library->is_rna, 1, 'is rna when extraction_type is DNA and transcript_strand is unstranded');

    $library->sample->extraction_type('total rna');
    $library->transcript_strand(undef);
    is($library->is_rna, 1, 'is rna when extraction_type is total rna and transcript_strand is UNDEF');

    $library->transcript_strand('unstranded');
    is($library->is_rna, 1, 'is rna when extraction_type is total rna and transcript_strand is unstranded');


};

done_testing();
