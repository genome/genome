#!/usr/bin/env genome-perl

#$ENV{UR_DBI_MONITOR_SQL}=1;

use above "Genome";

#use Test::More tests => 12;
use Test::More skip_all => 'many of the classes this tested have changed or no longer exist for this model';

my $model_id = 2661729970;

my $query_count;
Genome::DataSource::GMSchema->create_subscription(method => 'query',
                                                  callback => sub { $query_count++ },
                                                 );

my $load_count;
Genome::Model::Event->create_subscription(method => 'load',
                                          callback => sub { $load_count++},
                                      );

# By specifying a ref_seq_id, this will not include any events that need to join to
# the build table
my @events = Genome::Model::Command::InstrumentData::Assign->get(
                                                                 model_id => $model_id,
                                                                 instrument_data_id => 1
                                                             );
is(scalar(@events), 0 , "Genome::Model::Command::InstrumentData::Assign->get() correctly returns 0 items for model_id $model_id, chromosome 1");
is($query_count, 1, "get() generated 1 query");
ok($load_count, "at least one object was loaded by the get()");

$query_count = 0;
$load_count = 0;
@events = Genome::Model::Event::Build::ReferenceAlignment::MergeAlignments::Maq->get(model_id => $model_id, ref_seq_id => 1, event_status => 'Succeeded');
ok(scalar(@events), "Genome::Model::Event::Build::ReferenceAlignment::FindVariations->get() returned at least one event");
is($query_count, 0, "get() generated no queries");
is($load_count, 0, "and correctly loaded no objects");


$query_count = 0;
$load_count = 0;
# This one will include events that need to join to the build table
# Note that InstrumentData::Assign is now a real, logged event...
@events = Genome::Model::Command::InstrumentData::Assign->get(model_id => $model_id);
ok(scalar(@events) , "Genome::Model::Command::InstrumentData::Assign->get() correctly returned at least one event for model_id $model_id");
is($query_count, 1, "get() generated 1 query");
ok($load_count, "at least one object was loaded by the get()");
 
$query_count = 0;
$load_count = 0;
@events = Genome::Model::Event::Build::ReferenceAlignment::Solexa->get(model_id => $model_id);
ok(scalar(@events), "Genome::Model::Event::Build::ReferenceAlignment::Solexa->get() returned at least one event");
is($query_count, 0, "get() generated no queries");
is($load_count, 0, "and correctly loaded no objects");
