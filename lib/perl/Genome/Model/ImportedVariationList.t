#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
}
use above "Genome";
use Test::More tests => 6;

my $rsb = Genome::Model::Build::ImportedReferenceSequence->get(name => "NCBI-human-build36");
my $dbsnp = Genome::Model::ImportedVariationList->dbsnp_build_for_reference($rsb);
ok($dbsnp, "found dbsnp for reference " . $rsb->name);
ok($dbsnp->reference->is_compatible_with($rsb), "returned dbsnp build is compatible with reference");

$rsb = Genome::Model::Build::ImportedReferenceSequence->get(name => "g1k-human-build37");
$dbsnp = Genome::Model::ImportedVariationList->dbsnp_build_for_reference($rsb);
ok($dbsnp, "found dbsnp for reference " . $rsb->name);
ok($dbsnp->reference->is_compatible_with($rsb), "returned dbsnp build is compatible with reference");

$rsb = Genome::Model::Build::ImportedReferenceSequence->get(name => "GRCh37-lite-build37");
$dbsnp = Genome::Model::ImportedVariationList->dbsnp_build_for_reference($rsb);
ok($dbsnp, "found dbsnp for reference " . $rsb->name);
ok($dbsnp->reference->is_compatible_with($rsb), "returned dbsnp build is compatible with reference");

done_testing();
