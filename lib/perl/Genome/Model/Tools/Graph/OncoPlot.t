#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 2;
use File::Basename;
#use File::Spec my $dirName = dirname(__FILE__);
use File::Temp qw(tempfile);


# Creates temporary files
my $copyCatFileSample1 = new File::Temp( UNLINK => 0 );
my $copyCatFileSample2 = new File::Temp( UNLINK => 0 );
my $copyCatFileSample3 = new File::Temp( UNLINK => 0 );

my $variantFileSample1 = new File::Temp( UNLINK => 0 );
my $variantFileSample2 = new File::Temp( UNLINK => 0 );
my $variantFileSample3 = new File::Temp( UNLINK => 0 );

my $geneModelFile = new File::Temp( UNLINK => 0 );

my $inputFile = new File::Temp( UNLINK => 0 );


# Gene model file
my $geneFileContents = "chr\tstart\tend\tgene-name\n";
$geneFileContents .= "1\t2\t2000\tPine\n";
$geneFileContents .= "2\t5\t10\tLindel\n";
$geneFileContents .= "X\t5000\t10000000\tNewstead3\n";
#print $geneFileContents;
print $geneModelFile $geneFileContents;
close $geneModelFile;

# Sample variant files
my $sample1VariantData = "Gene\tChange\nPine\tSNP\n";
print $variantFileSample1 $sample1VariantData;
#print $sample1VariantData;
#print "\n";
close $variantFileSample1;

my $sample2VariantData = "Gene\tChange\nPine\tSNP\nLindel\tInsertion\nLindel\tSNP\n";
print $variantFileSample2 $sample2VariantData;
#print $sample2VariantData;
#print "\n";
close $variantFileSample2;

my $sample3VariantData = "Gene\tChange\nNewstead3\tFrameShift\n";
print $variantFileSample3 $sample3VariantData;
#print $sample3VariantData;
#print "\n";
close $sample3VariantData;

# Sample copycat data
my $copyCatSample1Data = "1\t3\t4\t10\t1\nX\t4\t5000\t10\t1\n";
print $copyCatFileSample1 $copyCatSample1Data;
close $copyCatFileSample1;

my $copyCatSample2Data = "2\t3\t4\t10\t2\n";
print $copyCatFileSample2 $copyCatSample2Data;
close $copyCatFileSample2;

my $copyCatSample3Data = "x\t5000\t10000000\t10\t1\nX\t5000\t10000000\t10\t3\n"; # won't match because the case is wrong
print $copyCatFileSample3 $copyCatSample3Data;
close $copyCatFileSample3;

# Inputs file : sample name; copycat file; variant file
my $inputs = join "\t", "SampleA", "$copyCatFileSample1", "$variantFileSample1";
$inputs .= "\n";
$inputs .= join "\t", "SampleB", "$copyCatFileSample2", "$variantFileSample2";
$inputs .= "\n";
$inputs .= join "\t", "SampleC", "$copyCatFileSample3", "$variantFileSample3\n";

#print $inputs;
print $inputFile $inputs;

my $outFile = new File::Temp( UNLINK => 0, SUFFIX => ".png");
my $genes = 'Pine,Lindel,Newstead3';
my $plot = Genome::Model::Tools::Graph::OncoPlot->create(
                                                         genes => $genes,
                                                         input => "$inputFile",
                                                         geneModel => "$geneModelFile",
                                                         outFile => "$outFile",
                                                         variantIndex => 2,
                                                         geneNameIndex => 1
                                                        );
ok($plot->execute,'Executed');
ok((-e "$outFile"), 'Image created');

#my $file = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Graph/heatmap-test-matrix.csv";
#print($file)                                                     