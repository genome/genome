#!/usr/bin/env genome-perl
#
#
# Tests all MC16s reports
#
#

use strict;
use warnings;

use above 'Genome';

use Genome::Model::Build::MetagenomicComposition16s::TestBuildFactory;
use Test::More;

use_ok('Genome::Model::Build::MetagenomicComposition16s::Reports') or die;

my ($build, $example_build) = Genome::Model::Build::MetagenomicComposition16s::TestBuildFactory->build_with_example_build_for_sanger;
my @amplicon_sets = $build->amplicon_sets;
my @example_amplicon_sets = $example_build->amplicon_sets;
ok(@amplicon_sets && @example_amplicon_sets, 'Got amplicon sets');
for ( my $i = 0; $i < @example_amplicon_sets; $i++ ) {
    for my $file_name (qw/ oriented_fasta_file oriented_qual_file classification_file /) {
        my $file = $example_amplicon_sets[$i]->$file_name;
        die "File ($file_name: $file) does not exist!" if not -s $file;
        Genome::Sys->create_symlink($file, $amplicon_sets[$i]->$file_name);
    }
}

# set some values
$build->amplicons_attempted(5);
is($build->amplicons_attempted, 5, 'amplicons attempted set');
$build->amplicons_processed(5);
is($build->amplicons_processed, 5, 'amplicons processed set');
$build->amplicons_processed_success(1);
is($build->amplicons_processed_success, 1, 'amplicons processed success set');

$build->amplicons_classified(4);
is($build->amplicons_classified, 4, 'amplicons classified set');
$build->amplicons_classified_success(4);
is($build->amplicons_classified_success, 4, 'amplicons classified success set');

# run
my $reports = Genome::Model::Build::MetagenomicComposition16s::Reports->create(input_build => $build);
ok($reports, 'create');
$reports->dump_status_messages(1);
ok($reports->execute, 'execute');

# verify
my @reports = glob($build->reports_directory.'/*');
is(@reports, 2, "Created 2 reports");
ok(-s $build->reports_directory.'/Summary/report.xml', 'Created summary report');
ok(-s $build->reports_directory.'/Summary/report.html', 'Created summary report html');
ok(-s $build->reports_directory.'/Composition/report.xml', 'Created composition report');

#print Data::Dumper::Dumper(map { $_->xml_string } $build->reports);<STDIN>;
done_testing();
