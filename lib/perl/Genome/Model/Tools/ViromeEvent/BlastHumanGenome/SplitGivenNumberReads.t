#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::ViromeEvent::BlastHumanGenome::SplitGivenNumberReads') or die;

#check testsuite files/dirs
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/Titanium17';
ok( -d $data_dir, "Testsuite data dir exists" ) or die;

my $run = 'Titanium17';
my $sample = $run.'_undecodable';
my $input_file = $sample.'.fa.cdhit_out.masked.goodSeq';
my $blast_dir = $sample.'.fa.cdhit_out.masked.goodSeq_HGblast';

#copy/link files
my $temp_dir = Genome::Sys->create_temp_directory();
ok( -d $temp_dir, "Created temp test dir" ) or die;

Genome::Sys->create_directory( $temp_dir."/$run" );
ok( -d $temp_dir."/$run", "Created temp run dir" ) or die;

Genome::Sys->create_directory( $temp_dir."/$run/$sample" );
ok( -d $temp_dir."/$run/$sample", "Created temp sample dir" ) or die;

symlink( $data_dir."/$sample/$input_file", $temp_dir."/$run/$sample/$input_file" );
ok( -l $temp_dir."/$run/$sample/$input_file", "Linked sample file" ) or die;

#create
my $c = Genome::Model::Tools::ViromeEvent::BlastHumanGenome::SplitGivenNumberReads->create(
    dir => $temp_dir."/$run/$sample",
    logfile => $temp_dir.'/log.txt',
    );

ok( $c, "Created blast human genome split-given-number-reads event" );

ok( $c->execute, "Successfully executed event" );

ok( -s $temp_dir."/$run/$sample/$blast_dir/$sample".'.fa.cdhit_out.masked.goodSeq_file0.fa', "Created pooled fasta file" );

done_testing();
