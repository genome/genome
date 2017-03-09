#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More tests => 5;

use above 'Genome';

use Genome::Test::Factory::Build;
use Genome::Test::Factory::Model::ReferenceAlignment;
use Genome::Utility::Test;

my $pkg = 'Genome::Model::Build::ReferenceAlignment::View::Coverage::Xml';
use_ok($pkg);

my $data_dir = Genome::Utility::Test->data_dir_ok($pkg, 'v1');
my $subject = setup_data($data_dir);

my $view_obj = Genome::Model::Build::ReferenceAlignment::View::Coverage::Xml->create(
    subject_id => $subject->id,
    aspects => ['status'],
);
isa_ok($view_obj,'Genome::Model::Build::ReferenceAlignment::View::Coverage::Xml');
my $xml = $view_obj->_generate_content;
ok($xml,'got xml content');


sub setup_data {
    my $data_dir = shift;
    my $m = Genome::Test::Factory::Model::ReferenceAlignment->setup_object();

    my $fl = Genome::FeatureList->__define__(
        name => 'Testing RefAlign Coverage View',
        content_type => 'exome',
        description => 'testing',
        format => 'true-BED',
        reference => $m->reference_sequence_build,
    );
    #$m->add_input( name => 'target_region_set_name', value_id => $fl->name, value_class_name => 'UR::Value::Text');
    #$m->add_input( name => 'region_of_interest_set_name', value_id => $fl->name, value_class_name => 'UR::Value::Text');

    my $subject = Genome::Test::Factory::Build->setup_object(model_id => $m->id);
    ok($subject, 'generated build subject');

    my $merged_result = Genome::InstrumentData::AlignmentResult::Merged->__define__(
        reference_build => $m->reference_sequence_build,
    );

    my $result = Genome::InstrumentData::AlignmentResult::Merged::CoverageStats->__define__(
        output_dir => File::Spec->join($data_dir, 'coverage_result'),
        wingspan_values => '0,500',
        minimum_depths => '1,5,10,15,20,30,40',
        region_of_interest_set => $fl,
        alignment_result => $merged_result,
    );
    $result->add_user( user => $subject, label => 'created');
    $result->add_user( user => $subject, label => 'uses');
    for my $type (qw(alignment alignment-v2 coverage)) {
        for my $wingspan (split ',', $result->wingspan_values) {
            my @depth = (($type eq 'coverage')? (1) : ());
            $result->add_metric(
                metric_name => join('_', join('-', $type, 'wingspan'), $wingspan, @depth, 'test-metric'),
                metric_value => 5,
            );
        }
    }

    return $subject;
}
