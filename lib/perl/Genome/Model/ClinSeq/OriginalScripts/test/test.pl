#!/usr/bin/env genome-perl
#Written by Malachi Griffith

#Load modules
use strict;
use warnings;
use above 'Genome'; #Makes sure I am using the local genome.pm and use libaries right above it (don't do this in a module)
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;


my $wgs_build_id = 114815422;
my $exome_build_id = 114752034;
my $rnaseq_build_id = 115909698;

my $wgs_build = Genome::Model::Build->get($wgs_build_id);
my $exome_build = Genome::Model::Build->get($exome_build_id);
my $rnaseq_build = Genome::Model::Build->get($rnaseq_build_id);

#my $wgs_common_name = $wgs_build->subject->patient->common_name;
#my $exome_common_name = $exome_build->subject->patient->common_name;
#my $rnaseq_common_name = $rnaseq_build->subject->patient->common_name;

my $wgs_common_name = '';
my $exome_common_name = '';
my $rnaseq_common_name = '';


my $wgs_patient_name = $wgs_build->subject->patient->name;
my $exome_patient_name = $exome_build->subject->patient->name;
my $rnaseq_patient_name = $rnaseq_build->subject->patient->name;

my $wgs_subject_name = $wgs_build->subject->name;
my $exome_subject_name = $exome_build->subject->name;
my $rnaseq_subject_name = $rnaseq_build->subject->name;

print BLUE, "\n\nwgs_common_name = $wgs_common_name\twgs_subject_name = $wgs_subject_name\twgs_patient_name = $wgs_patient_name\nexome_common_name = $exome_common_name\texome_subject_name = $exome_subject_name\texome_patient_name = $exome_patient_name\nrnaseq_common_name = $rnaseq_common_name\trnaseq_subject_name = $rnaseq_subject_name\trnaseq_patient_name = $rnaseq_patient_name\n\n", RESET;

my @names = ($wgs_common_name, $exome_common_name, $rnaseq_common_name, $wgs_patient_name, $exome_patient_name, $rnaseq_patient_name, $wgs_subject_name, $exome_subject_name, $rnaseq_subject_name);

my $final_name;

foreach my $name (@names){
  if ($name){
    $final_name = $name;
    last();
  }
}

print BLUE, "\n\nFinal name = $final_name\n\n", RESET;

exit();

