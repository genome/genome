
use strict;
use warnings;
use above 'Genome';

use Test::More tests => 6;

BEGIN
{
    use_ok("Genome::Model::ImportedAnnotation");
}

my $model = Genome::Model::ImportedAnnotation->get(2771411739);
isa_ok($model, "Genome::Model::ImportedAnnotation");
ok($model->name, 'human.imported-annotation-NCBI-human-36');
my $return =  eval {my $build = $model->build_by_version(0);};
ok(!$return, 'Fetching build by version 0 fails');
my $build = $model->build_by_version('54_36p_v2');
ok($build, "Fetching build by version '54_36p_v2'");
isa_ok($build, "Genome::Model::Build::ImportedAnnotation");
