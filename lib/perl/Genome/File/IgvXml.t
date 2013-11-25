#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 10;

#test class creation
my $path = Genome::Sys->create_temp_file_path;
my $igv_file = Genome::File::IgvXml->create(id => $path, font_size => 10, genome => 'b37');
ok($igv_file, "created a IgvXml object");
is($igv_file->path,$path,"path assigned correctly");
is($igv_file->font_size,10, "class variable set ok on create");

#tests on functions
my $bam = "/some/fake.bam";
$igv_file->add_resource($bam);

is("    <Resource path=\"$bam\"/>", $igv_file->_resources->[0], "resource added correctly");
$igv_file->add_data_track($bam);
is($bam,$igv_file->_data_tracks->[0], "track added correctly");
$igv_file->add_feature_track($bam);
is($bam, $igv_file->_feature_tracks->[0], "feature track added correctly");

my $expected_panel1 = <<XML;
  <Panel name="Panel 1">
$bam
  </Panel>
XML

is($expected_panel1, $igv_file->make_panel($bam), "panel created with autonaming");

my $expected_panel2 = <<XML;
  <Panel name="Test Panel">
$bam
  </Panel>
XML

is($expected_panel2, $igv_file->make_panel($bam, "Test Panel"), "panel created with custom name");

my $path2 = Genome::Sys->create_temp_file_path;
my $igv_file2 = Genome::File::IgvXml->create(id => $path2, font_size => 10, genome => 'b37');
my $fake_bam = "/some/fake.bam";
my $fake_bed = "/some/fake.bed";
my $fake_junctions = "/some/junctions.txt"; #no idea what a junction file is

$igv_file2->add_bam_track(file => $fake_bam);
$igv_file2->add_bed_track(file => $fake_bed);
$igv_file2->add_junction_track(file => $fake_junctions);

my $expected_xml = <<XML;
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="b37" locus="12:25398182-25398361" version="4">
  <Resources>
    <Resource path="https://gscweb.gsc.wustl.edu/some/fake.bam"/>
    <Resource path="https://gscweb.gsc.wustl.edu/some/fake.bed"/>
    <Resource path="https://gscweb.gsc.wustl.edu/some/junctions.txt"/>
  </Resources>
  <Panel name="Panel 1">
    <Track altColor="0,0,178" autoScale="true" color="175,175,175" colorScale="ContinuousColorScale;0.0;9062.0;255,255,255;175,175,175" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="https://gscweb.gsc.wustl.edu/some/fake.bam_coverage" name="/some/fake.bam Coverage" showDataRange="true" visible="true">
      <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="100" minimum="0" type="LINEAR"/>
    </Track>
    <Track altColor="0,0,178" color="0,0,178" colorOption="UNEXPECTED_PAIR" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" id="https://gscweb.gsc.wustl.edu/some/fake.bam" name="/some/fake.bam Reads" showDataRange="true" sortByTag="" visible="true"/>
  </Panel>
  <Panel name="FeaturePanel">
    <Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" showDataRange="true" sortable="false" visible="true"/>
    <Track altColor="0,0,178" color="0,0,178" colorScale="ContinuousColorScale;0.0;160.0;255,255,255;0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="35" id="b37 Genes" name="Gene" renderer="BASIC_FEATURE" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>
    <Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="45" id="https://gscweb.gsc.wustl.edu/some/fake.bed" name="/some/fake.bed" renderer="BASIC_FEATURE" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>    <Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="60" id="https://gscweb.gsc.wustl.edu/some/junctions.txt" name="/some/junctions.txt" showDataRange="true" visible="true" windowFunction="count"/>
  </Panel>
  <PanelLayout dividerFractions="0.000,0.450"/>
</Session>
XML

is($igv_file2->xml, $expected_xml, "XML with Bam file, junction file and bed file as expected.");

my $path3 = Genome::Sys->create_temp_file_path;
my $igv_file3 = Genome::File::IgvXml->create(id => $path3, font_size => 10, genome => 'b37');
$igv_file3->add_bam_track(file => $fake_bam);
$igv_file3->add_bed_track(file => $fake_bed);
$igv_file3->add_junction_track(file => $fake_junctions);

my $expected_xml2 = <<XML;
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="b37" locus="1:2-5" version="4">
  <Resources>
    <Resource path="https://gscweb.gsc.wustl.edu/some/fake.bam"/>
    <Resource path="https://gscweb.gsc.wustl.edu/some/fake.bed"/>
    <Resource path="https://gscweb.gsc.wustl.edu/some/junctions.txt"/>
  </Resources>
  <Panel name="Panel 1">
    <Track altColor="0,0,178" autoScale="true" color="175,175,175" colorScale="ContinuousColorScale;0.0;9062.0;255,255,255;175,175,175" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="https://gscweb.gsc.wustl.edu/some/fake.bam_coverage" name="/some/fake.bam Coverage" showDataRange="true" visible="true">
      <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="100" minimum="0" type="LINEAR"/>
    </Track>
    <Track altColor="0,0,178" color="0,0,178" colorOption="UNEXPECTED_PAIR" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" id="https://gscweb.gsc.wustl.edu/some/fake.bam" name="/some/fake.bam Reads" showDataRange="true" sortByTag="" visible="true"/>
  </Panel>
  <Panel name="FeaturePanel">
    <Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" showDataRange="true" sortable="false" visible="true"/>
    <Track altColor="0,0,178" color="0,0,178" colorScale="ContinuousColorScale;0.0;160.0;255,255,255;0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="35" id="b37 Genes" name="Gene" renderer="BASIC_FEATURE" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>
    <Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="45" id="https://gscweb.gsc.wustl.edu/some/fake.bed" name="/some/fake.bed" renderer="BASIC_FEATURE" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>    <Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="60" id="https://gscweb.gsc.wustl.edu/some/junctions.txt" name="/some/junctions.txt" showDataRange="true" visible="true" windowFunction="count"/>
  </Panel>
  <PanelLayout dividerFractions="0.000,0.450"/>
</Session>
XML

is($igv_file3->xml("1:2-5"), $expected_xml2,"XML with Bam file, junction file and bed file with custom locus as expected.");
