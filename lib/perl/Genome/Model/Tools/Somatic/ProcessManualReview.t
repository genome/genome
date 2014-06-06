#!/usr/bin/env genome-perl

use strict;
use warnings;
use above "Genome";
use Test::More tests => 4;
use File::Spec;

my $variantMatcherOut = new File::Temp( UNLINK => 0 );

my $variantFile = File::Spec->catfile(
	$ENV{GENOME_TEST_INPUTS}, "Genome-Model",
	"ProcessManualReview",    "Snvs.indels.annotated"
);
my $reviewedFile = File::Spec->catfile(
	$ENV{GENOME_TEST_INPUTS}, "Genome-Model",
	"ProcessManualReview",    "Reviewed.bed"
);
use_ok('Genome::Model::Tools::Somatic::ProcessManualReview');
#my $analysis = Genome::Model::Tools::Somatic::ProcessManualReview->create();

my $analysis = Genome::Model::Tools::Somatic::ProcessManualReview->create(
	statuses      => "O,S",
	reviewed_file => $reviewedFile,
	variant_file  => $variantFile,
	output_file   => "$variantMatcherOut"
);

#my $plot = Genome::Model::Tools::Graph::OncoPlot->create(
#                                                         genes => $genes,
#                                                         input => "$inputFile",
#                                                         geneModel => "$geneModelFile",
#                                                         outFile => "$outFile",
#                                                         variantIndex => 2,
#                                                         geneNameIndex => 1
#                                                        );
 

ok( $analysis->execute, "Execution successful" );
ok( ( -e "$variantMatcherOut" ), "Ouput file created" );

# Test the output
my $line_count = 0;
my $word_count = 0;

while ( my $line = <$variantMatcherOut> ) {
	$line_count++;
	my @words_on_this_line = split /\t/, $line;
	$word_count += scalar(@words_on_this_line);
}
close $variantMatcherOut;

ok( ( $line_count == 11 && $word_count == 374 ), "Output looks great!" );
