#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

my $lims_class = 'Genome::Site::TGI::Synchronize::Classes::InstrumentDataAnalysisProjectBridge';
use_ok($lims_class) or die;

my $entity_name = $lims_class->entity_name;
is($entity_name, 'analysis project instrument data', 'entity name');
is($lims_class->genome_class_for_comparison, 'Genome::Site::TGI::Synchronize::Classes::AnalysisProjectInstrumentData', 'genome class for comparison');
is($lims_class->genome_class_for_create, 'Genome::Config::AnalysisProject::InstrumentDataBridge', 'genome class for create');

my @properties_to_copy = $lims_class->properties_to_copy;
ok(@properties_to_copy, 'properties to copy');
my $i = -10;
my %properties = ( analysis_project_id => --$i, instrument_data_id => --$i );
my $lims_object = $lims_class->__define__(%properties);
ok($lims_object, "define lims $entity_name object");

my $genome_object = $lims_object->create_in_genome;
ok($genome_object, "create genome $entity_name object");
isa_ok($genome_object, 'Genome::Config::AnalysisProject::InstrumentDataBridge');

is($genome_object->analysis_project_id, $properties{analysis_project_id}, 'analysis_project_id matches');
is($genome_object->instrument_data_id, $properties{instrument_data_id}, 'instrument_data_id matches');

done_testing();
