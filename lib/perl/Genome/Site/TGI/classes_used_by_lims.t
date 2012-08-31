#!/usr/bin/env genome-perl 

# This tests Genome classes used in LIMS (GSC). They should be tested using the LIMS verison of perl.
# There is a reciprocal test in LIMS for GSC classes used by Genome. See ebelter or idas for more info.

use warnings;
use strict;

use Test::More;

my @apipe_classes_used_by_lims = (qw/
    Genome::InstrumentData::Report
    Genome::InstrumentData::Solexa
    Genome::Library
    Genome::Model
    Genome::Model::Build
    Genome::Model::Build::ImportedVariationList
    Genome::Model::Build::ReferenceAlignment
    Genome::Model::Build::ReferenceSequence
    Genome::Model::Input
    Genome::Model::Metric
    Genome::Model::ReferenceAlignment::Command::CreateMetrics::CompareSnps
    Genome::Model::SomaticVariation
    Genome::Model::Tools::ContaminationScreen::3730
    Genome::Model::Tools::Fastqc::GenerateReports
    Genome::Model::Tools::Joinx::SnvConcordanceByQuality
    Genome::Model::Tools::Picard
    Genome::Model::Tools::Picard::CollectGcBiasMetrics
    Genome::Model::Tools::Picard::EstimateLibraryComplexity
    Genome::Model::Tools::Picard::FastqToSam
    Genome::Model::Tools::Picard::ReplaceSamHeader
    Genome::Model::Tools::Picard::SamToFastq
    Genome::Model::Tools::Sam::SortBam
    Genome::ModelGroup
    Genome::SoftwareResult
    /);

diag("Perl Version: $]");
for my $class ( @apipe_classes_used_by_lims ) {
    my $cmd = "/gsc/bin/perl -MGenome -e '$class->class;'";
    #print $cmd, "\n";
    system "$cmd";
    is($?, 0, "$class is usable by LIMS");
}

done_testing( scalar(@apipe_classes_used_by_lims) );
exit;

