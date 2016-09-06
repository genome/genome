#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above "Genome";

use Test::More tests => 2;

subtest 'use and roles' => sub{
    plan tests => 2;

    use_ok('Genome::Library') or die;
    ok(Genome::Library->does('Genome::Role::Notable'), 'Library does Genome::Role::Notable');

};

my $library;
subtest 'create' => sub{
    plan tests => 5;

    my $sample = Genome::Sample->create(name => '__TEST_SAMPLE__');
    ok($sample, 'create sample');

    $library = Genome::Library->create(
        sample => $sample,
        name => $sample->name . "-extlibs",
        original_insert_size => '1kb',
        library_insert_size => '300-500',
        protocol => 'karate chop',
        transcript_strand => 'unstranded',
    );
    isa_ok($library, 'Genome::Library');
    is($library->name, $sample->name . "-extlibs", "name is what is expected");
    ok($library->is_rna, 'is_rna true when transcript strand is set');

    my $commit = eval{ UR::Context->commit; };
    ok($commit, 'commit');

};

done_testing();
