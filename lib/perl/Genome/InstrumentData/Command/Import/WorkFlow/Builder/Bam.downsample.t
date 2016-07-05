#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More tests => 6;

my $class = 'Genome::InstrumentData::Command::Import::WorkFlow::Builder::Bam';
use_ok($class) or die;

my @source_files = (qw/ in.bam /);
my $analysis_project = Genome::Config::AnalysisProject->__define__(name => 'TEST-AnP');
ok($analysis_project, 'define analysis project');
my $library = Genome::Library->__define__(name => 'TEST-sample-libs', sample => Genome::Sample->__define__(name => 'TEST-sample'));
ok($library, 'define library');

my $builder = Genome::InstrumentData::Command::Import::WorkFlow::Builder->create(
    work_flow_inputs => Genome::InstrumentData::Command::Import::Inputs::Factory->create(
        analysis_project => $analysis_project,
    )->from_params({
            source_paths => \@source_files,
            entity_params => { instdata => { downsample_ratio => .25, }, }, 
        }),
);
isa_ok($builder, $class);
my $wf = $builder->build_workflow;
ok($wf, 'build_workflow');
my $expected_xml_file = __FILE__;
$expected_xml_file =~ s/t$/xml/;
my $expected_xml = Genome::Sys->read_file($expected_xml_file);
is($wf->get_xml, $expected_xml, 'xml matches');

done_testing();
