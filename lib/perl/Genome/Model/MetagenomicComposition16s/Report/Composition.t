#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Model::MetagenomicComposition16s::Test;
use Test::More;

use_ok('Genome::Model::MetagenomicComposition16s::Report::Composition') or die;

my $model = Genome::Model::MetagenomicComposition16s::Test->model_for_sanger;
ok($model, 'got mc16s sanger model');
my $example_build = Genome::Model::MetagenomicComposition16s::Test->example_build_for_model($model);
ok($example_build, 'got example build');

my $generator = Genome::Model::MetagenomicComposition16s::Report::Composition->create(build_id => $example_build->id);
ok($generator, 'create');
my $report = $generator->generate_report;
ok($report, 'generate report');

#print $report->xml_string."\n";
#print $example_build->classification_file_for_set_name('');<STDIN>;
done_testing();
exit;

