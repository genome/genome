#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above "Genome";

use Test::More;

use_ok('Genome::Library') or die;

my $sample = Genome::Sample->create(name => '__TEST_SAMPLE__');
ok($sample, 'create sample');
my $library = Genome::Library->create(
    sample => $sample,
    name => $sample->name . "-extlibs",
    original_insert_size => '1kb',
    library_insert_size => '300-500',
    protocol => 'karate chop',
    transcript_strand => 'unstranded',
);
ok($library, 'create library');
isa_ok($library, 'Genome::Library');
ok($library->does('Genome::Role::Notable'), 'Library does Genome::Role::Notable');
is($library->name, $sample->name . "-extlibs", "name is what is expected");

ok($library->is_rna, 'is_rna true when transcript strand is set');

my $commit = eval{ UR::Context->commit; };
ok($commit, 'commit');

done_testing();
