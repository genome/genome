#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::Exception;
use Test::More;

my $pkg = 'Genome::VariantReporting::Suite::Vep::SpliceNonsynonymousList';
use_ok($pkg);

subtest 'is_splice_site' => sub {
    ok(
        Genome::VariantReporting::Suite::Vep::SpliceNonsynonymousList::is_splice_site('splice_acceptor_variant'), 
        "('splice_acceptor_variant') is splice site"
    );
    ok(
        Genome::VariantReporting::Suite::Vep::SpliceNonsynonymousList::is_splice_site('splice_acceptor_variant', 'not_splice_site'), 
        "('splice_acceptor_variant', 'not_splice_site') is splice site"
    );
    ok(
        !Genome::VariantReporting::Suite::Vep::SpliceNonsynonymousList::is_splice_site('not_splice_site'), 
        "('not_splice_site') is not splice site"
    );
};

subtest 'is_non_synonymous' => sub {
    ok(
        Genome::VariantReporting::Suite::Vep::SpliceNonsynonymousList::is_non_synonymous('transcript_ablation'), 
        "('transcript_ablation') is non synonymous"
    );
    ok(
        Genome::VariantReporting::Suite::Vep::SpliceNonsynonymousList::is_non_synonymous('transcript_ablation', 'not_non_synonymous'), 
        "('transcript_ablation', 'not_non_synonymous') is non synonymous"
    );
    ok(
        !Genome::VariantReporting::Suite::Vep::SpliceNonsynonymousList::is_non_synonymous('not_non_synonymous'), 
        "('not_non_synonymous') is not non synonymous");
};

done_testing;
