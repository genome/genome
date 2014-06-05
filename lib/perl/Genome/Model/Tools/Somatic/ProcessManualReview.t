#!/usr/bin/env genome-perl

use strict;
use warnings;
use Test::More tests => 3;
#use ProcessManualReview;
use File::Spec;
use above 'Genome';
my $variantMatcherOut = new File::Temp( UNLINK => 0 );

my $variantFile = File::Spec->catfile(
	$ENV{GENOME_TEST_INPUTS}, "Genome-Model",
	"ProcessManualReview",    "Snvs.indels.annotated"
);
my $reviewedFile = File::Spec->catfile(
	$ENV{GENOME_TEST_INPUTS}, "Genome-Model",
	"ProcessManualReview",    "Reviewed.bed"
);
my $analysis = Genome::Model::Tools::Somatic::ProcessManualReview->create(
	statuses      => "O,S",
	reviewed_file => $reviewedFile,
	variant_file  => $variantFile,
	output_file   => "$variantMatcherOut",
);

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
