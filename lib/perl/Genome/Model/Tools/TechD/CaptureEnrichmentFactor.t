#!/usr/bin/env genome-perl

# Simple test of CaptureEnrichmentFactor class.

use strict;
use warnings;
use above "Genome";

use Test::More tests => 5;

use_ok('Genome::Model::Tools::TechD::CaptureEnrichmentFactor') or die "test cannot continue...";


my $myEF = Genome::Model::Tools::TechD::CaptureEnrichmentFactor->execute(
                                        capture_unique_bp_on_target    =>  6_795_966_000,
                                        capture_duplicate_bp_on_target =>    834_157_200,
                                        capture_total_bp               => 13_613_105_800,
                                        target_total_bp                =>     45_203_256,
                                        genome_total_bp                =>  2_858_012_809,
                                       );

ok($myEF->result, "CaptureEnrichmentFactor executed.") or die "test cannot continue...";

my $theoretical_max_enrichment_factor  = $myEF->theoretical_max_enrichment_factor();
ok($theoretical_max_enrichment_factor eq 63.2) or die "theoretical_max_enrichment_factor should equal 6.2, but $theoretical_max_enrichment_factor is returned";

my $unique_on_target_enrichment_factor = $myEF->unique_on_target_enrichment_factor();
ok($unique_on_target_enrichment_factor eq 31.6) or die "unique_on_target_enrichment_factor should equal 31.6, but $unique_on_target_enrichment_factor is returned";

my $total_on_target_enrichment_factor  = $myEF->total_on_target_enrichment_factor();
ok($total_on_target_enrichment_factor eq 35.4) or die "total_on_target_enrichment_factor should equal 35.4, but $total_on_target_enrichment_factor is returned";

__END__
