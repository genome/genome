#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Model::MetagenomicComposition16s::Test;
use Test::More;

use_ok('Genome::Model::MetagenomicComposition16s::Report::Summary') or die;

my $model = Genome::Model::MetagenomicComposition16s::Test->model_for_sanger;
ok($model, 'Got model') or die;
my $build = Genome::Model::MetagenomicComposition16s::Test->example_build_for_model($model);
ok($build, 'Got build') or die;

# check that it works w/o metrics
my $generator = Genome::Model::MetagenomicComposition16s::Report::Summary->create(build_id => $build->id);
ok($generator, 'Created generator');
my $report = $generator->generate_report;
ok($report, 'Generated report');
#print $report->xml_string."\n";

# set metrics
$build->amplicons_attempted(5);
$build->amplicons_processed(4);
$build->amplicons_processed_success('0.80');
$build->amplicons_classified(4);
$build->amplicons_classified_success(1);
$build->reads_attempted(30);
$build->reads_processed(17);
$build->reads_processed_success('0.57');

# real report
$generator = Genome::Model::MetagenomicComposition16s::Report::Summary->create(build_id => $build->id);
ok($generator, 'Created generator');
$report = $generator->generate_report;
ok($report, 'Generated report');

#print $report->xml_string."\n";
#print $build->classification_file_for_set_name('');<STDIN>;
done_testing();
exit;

