use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above "Genome";

use Genome::Test::Factory::Model::ReferenceAlignment;
use Genome::Test::Factory::InstrumentData::Solexa;

use constant NUM_INSTRUMENT_DATA => 5;
use Test::More tests => 8 + NUM_INSTRUMENT_DATA;

use_ok('Genome::Model::Command::SplitOnLibrary');

my $ra_model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object();
isa_ok($ra_model, 'Genome::Model', 'generated model');

my @l = map {
    Genome::Test::Factory::Library->setup_object(
        sample_id => $ra_model->subject_id
    )
} 1..3;
ok(scalar(@l), 'generated libraries');

my @i = map {
    Genome::Test::Factory::InstrumentData::Solexa->setup_object(
        library_id => $l[$_ % scalar(@l)]
    )
} 1..NUM_INSTRUMENT_DATA;
ok(scalar(@i), 'generated instrument data');


map { $ra_model->add_instrument_data($_) } @i;

my $cmd = Genome::Model::Command::SplitOnLibrary->create(model => $ra_model);
isa_ok($cmd, 'Genome::Model::Command::SplitOnLibrary', 'generated command');

ok($cmd->execute, 'executed command');

my @new_models = $cmd->new_models;
is(scalar(@new_models), scalar(@l), 'generated one model per library');

my %seen_library_ids;
for my $m (@new_models) {
    my @model_instrument_data = $m->instrument_data;
    my $library_id =  $model_instrument_data[0]->library_id;
    $seen_library_ids{$library_id}++;

    for my $i (@model_instrument_data) {
        is($library_id, $i->library_id, 'library id matches for all instrument data in model');
    }
}

is(scalar(keys %seen_library_ids), scalar(@l), 'made one model per library');


