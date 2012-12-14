#!/usr/bin/env genome-perl
package Genome::Model::Tools::Pcap::Phd::Test;
use above 'Genome';
use Genome::Model::Tools::Pcap::Phd;
use Test::More tests => 1;

my $po = Genome::Model::Tools::Pcap::Phd->new(input_directory => $ENV{GENOME_TEST_INPUTS} . '/Genome-Assembly-Pcap/phd_dir');

my $phd = $po->get_phd("L25990P6007H3.g1.phd.1");
is($phd->name,"L25990P6007H3.g1","Phd name survives creation");
