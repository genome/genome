#!/usr/bin/env genome-perl

use above "Genome";
use Test::More skip_all => 'needs to have more reasonably sized test data before this runs in the suite, yell at adukes if you see this message a lot';


my $model = Genome::Model->get(name => "HMPZ-764224817-700024417.HMP_mWGS_100702_v1.0");
my $dd = $model->last_succeeded_build->data_directory();
system "ls $dd";

my @bams = ("$dd/metagenomic_alignment1.bam", "$dd/metagenomic_alignment2.bam");

for (@bams){
    ok(-e $_, "bam to merge exists");
}

$output_file = "/tmp/merged_bam.bam";

SKIP{
print join("\n", @INC)."\n";
ok(Genome::Model::Tools::Sam::SortAndMergeSplitReferenceAlignments->execute(
        input_files => \@bams,
        output_file => $output_file,
    )->result, 'executed sort and merge bam');

ok (-s $output_file, "output file has size");
}
