#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;
use Cwd;

#check test suite dir/files
my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-AssemblReads-Pcap/Proteus_penneri_ATCC_35198-1.0_080509.pcap';
ok (-d $test_dir, "Test suite directory exists");

#create temp test dir
my $temp_dir = Genome::Sys->create_temp_directory();
ok (-d $temp_dir, "Created temp test directory");

#parts of pcap must run in assembly directory
my $pre_test_dir = cwd();

#create/execute tool
my $obj = Genome::Model::Tools::Pcap::Assemble->create (
     project_name       => 'Proteus_penneri_ATCC_35198',
     disk_location      => $temp_dir,#$ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-AssemblReads-Pcap',
     parameter_setting  => 'RELAXED',
     assembly_version   => '1.0',
     assembly_date      => '080509',
     read_prefixes      => 'PPBA',
     pcap_run_type      => 'NORMAL',
    );
ok ($obj, "Instance of pcap-assemble tool created"); 
ok ($obj->_project_path(), "Object set project path");
ok ($obj->create_project_directories, "Created project dirs");

#check is assembly directory created
ok (-d $temp_dir.'/Proteus_penneri_ATCC_35198-1.0_080509.pcap', "Temp assembly dir created");
#link/copy input files into assembly directory
foreach (qw/ PPBA.fasta.gz PPBA.fasta.qual.gz /) {
    ok (-s $test_dir."/input/$_", "Test suite $_ input file exists");
    symlink ($test_dir."/input/$_", $temp_dir."/Proteus_penneri_ATCC_35198-1.0_080509.pcap/input/$_");
    ok (-l $temp_dir."/Proteus_penneri_ATCC_35198-1.0_080509.pcap/input/$_", "Input $_ file linked");
}

#test methods .. really should just run execute but tool runs off config file
#TODO - have test write config file and run execute with it
ok($obj->copy_test_data_set, "test data set copied");
ok($obj->create_pcap_input_fasta_fof, "pcap fasta fof created");
ok($obj->create_constraint_file, "constraint file created successfully");
ok($obj->resolve_pcap_run_type, "pcap run type resolved");
ok($obj->run_pcap, "test pcap.rep");
ok($obj->run_bdocs, "test bdocs.rep");
ok($obj->run_bclean, "test bclean.rep");
ok($obj->run_bcontig, "test bcontig.rep");
ok($obj->check_for_results_file, "test check for results file");
ok($obj->run_bconsen, "test bconsen.test");
ok($obj->run_bform, "test bform.rep");
ok($obj->create_gap_file, "test create_gap_file");
ok($obj->create_agp_file, "test create_agp_file");
ok($obj->create_sctg_fa_file, "test create_sctg_fa_file");
ok($obj->add_wa_tags_to_ace, "test add WA tags to ace");

#retrun to original dir to clean up test dir before existing
chdir $pre_test_dir;

done_testing();
