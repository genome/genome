use strict;
use warnings;

use Genome::Test::Factory::Sample qw();
use Genome::Test::Factory::Model::ReferenceAlignment qw();

use Test::More tests => 6;

my $common_name_prop = Genome::Sample->__meta__->properties(property_name => 'common_name');
is($common_name_prop->via, 'attributes', 'common_name is via attributes');

my $models_prop = Genome::Sample->__meta__->properties(property_name => 'models');
is($models_prop->reverse_as, 'subject', 'models is reverse_as subject');

my $sample = Genome::Test::Factory::Sample->setup_object(
    common_name => 'foo',
);
my $model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(
    subject_id => $sample->id,
);

is_deeply([$sample->models], [$model], 'associated model with sample');

my $copy = $sample->copy();

is($copy->common_name, 'foo', 'copy has same common_name');
isnt($sample->id, $copy->id, 'copy has a different id');
is(scalar(() = $copy->models), 0, 'copy has no models associted with it');
