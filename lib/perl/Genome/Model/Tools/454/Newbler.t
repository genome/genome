#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More skip_all => 'times out in Jenkins and causes failures, needs to be refactored to not take so long';
use File::Temp;
use File::Path;
use Data::Dumper;

BEGIN {
    my $archos = `uname -a`;
    if ($archos !~ /64/) {
        plan skip_all => "Must run from 64-bit machine";
    }
    plan tests => 66;

    use_ok( 'Genome::Model::Tools::454::Newbler::NewMapping');
    use_ok( 'Genome::Model::Tools::454::Newbler::NewAssembly');
    use_ok( 'Genome::Model::Tools::454::Newbler::SetRef');
    use_ok( 'Genome::Model::Tools::454::Newbler::AddRun');
    use_ok( 'Genome::Model::Tools::454::Newbler::RemoveRun');
    use_ok( 'Genome::Model::Tools::454::Newbler::RunProject');
    use_ok( 'Genome::Model::Tools::454::Newbler::RunMapping');
    use_ok( 'Genome::Model::Tools::454::Newbler::RunAssembly');
};

my $data_dir = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-454-Newbler';
#my $expected_path = Genome::Config::get('sw') . '/454/newbler/applicationsBin/';
#my $expected_path = Genome::Config::get('sw') . '/454/offInstrumentApps-2.0.00.20-64/bin';

my @version_subdirs = qw/ offInstrumentApps mapasm454_source /;

foreach my $subdir (@version_subdirs) {

    my ($version, $expected_path);
    if ($subdir =~ /offInstrumentApps/) {
	$version = '2.0.00.20-1';
	$expected_path = Genome::Config::get('sw') . '/454/'.$subdir.'-'.$version.'/bin';
    }
    if ($subdir =~ /mapasm454_source/) {
	$version = '10282008';
	$expected_path = Genome::Config::get('sw') . '/454/'.$subdir.'_'.$version.'/applicationsBin';
    }

    my $ref_seq_dir = $data_dir . '/refseq-for-test';
    my @fasta_files = glob($ref_seq_dir .'/11.fasta');
    is(scalar(@fasta_files),1,'correct number of test fasta files found');

    opendir(DATA,$data_dir) || die ("Can not open data directory '$data_dir'");
    my @runs = grep { !/^\./  } readdir DATA;
    closedir DATA;
    is(scalar(@runs),1,'correct number of test runs found');
    my $run_name = $runs[0];

    my $run_dir = $data_dir .'/'. $run_name;
    opendir(RUN,$run_dir) || die ("Can not open run directory '$run_dir'");
    my @sff_files = map { $run_dir.'/'.$_ } grep { /\.sff$/  } readdir RUN;
    closedir RUN;
    is(scalar(@sff_files),4,'correct number of test sff files found');

    my $mapping_project_dir = File::Temp::tempdir(CLEANUP => 1);
    ok(rmtree($mapping_project_dir),"removed tmp directory: $mapping_project_dir");
    my $new_mapping = Genome::Model::Tools::454::Newbler::NewMapping->create(
									     dir => $mapping_project_dir,
									     version => $version,
									     version_subdirectory => $subdir,
									     );
    isa_ok($new_mapping,'Genome::Model::Tools::454::Newbler::NewMapping');
    ok($new_mapping->execute,'execute newbler newMapping');
    ok(-d $mapping_project_dir,'newbler mapping directory created');
    my $set_ref = Genome::Model::Tools::454::Newbler::SetRef->create(
								     dir => $mapping_project_dir,
								     reference_fasta_files => \@fasta_files,
								     version => $version,
								     version_subdirectory => $subdir,
								     );
    isa_ok($set_ref,'Genome::Model::Tools::454::Newbler::SetRef');
    ok($set_ref->execute,'execute newbler setRef');
    
    my $assembly_project_dir = File::Temp::tempdir(CLEANUP => 1);
    ok(rmtree($assembly_project_dir),"removed tmp directory: $assembly_project_dir");
    my $new_assembly = Genome::Model::Tools::454::Newbler::NewAssembly->create(
									       dir => $assembly_project_dir,
									       version => $version,
									       version_subdirectory => $subdir,
									       );
    isa_ok($new_assembly,'Genome::Model::Tools::454::Newbler::NewAssembly');
    ok($new_assembly->execute,'execute newbler newAssembly');
    ok(-d $assembly_project_dir,'newbler assembly directory created');

    my @dirs = ($mapping_project_dir, $assembly_project_dir);
    foreach my $dir (@dirs) {
	my $add_run = Genome::Model::Tools::454::Newbler::AddRun->create(
									 dir => $dir,
									 runs => \@sff_files,
									 version => $version,
									 version_subdirectory => $subdir,
									 );
	isa_ok($add_run,'Genome::Model::Tools::454::Newbler::AddRun');
	ok($add_run->execute,'execute newbler addRun');
	
	my $run_project = Genome::Model::Tools::454::Newbler::RunProject->create(
										 dir => $dir,
										 version => $version,
										 version_subdirectory => $subdir,
										 );
	isa_ok($run_project,'Genome::Model::Tools::454::Newbler::RunProject');
	ok($run_project->execute,'execute newbler runProject');
    }

    my $mapping_dir = File::Temp::tempdir(CLEANUP => 1);
    ok(rmtree($mapping_dir),"removed tmp directory: $mapping_dir");
    my $run_mapping = Genome::Model::Tools::454::Newbler::RunMapping->create(
									     mapping_dir => $mapping_dir,
									     sff_files => \@sff_files,
									     ref_seq => $fasta_files[0],
									     version => $version,
									     version_subdirectory => $subdir,
									     );
    isa_ok($run_mapping,'Genome::Model::Tools::454::Newbler::RunMapping');
    ok($run_mapping->execute,'execute newbler runMapping');

    my $assembly_dir = File::Temp::tempdir(CLEANUP => 1);
    ok(rmtree($assembly_dir),"removed tmp directory: $assembly_dir");
    my $run_assembly = Genome::Model::Tools::454::Newbler::RunAssembly->create(
									       assembly_dir => $assembly_dir,
									       sff_files => \@sff_files,
									       version => $version,
									       version_subdirectory => $subdir,
									       );
    isa_ok($run_assembly,'Genome::Model::Tools::454::Newbler::RunAssembly');
    ok($run_assembly->execute,'execute newbler runAssembly');
    
    my $version_run_assembly = Genome::Model::Tools::454::Newbler::RunAssembly->create (
											assembly_dir => $assembly_dir,
											sff_files => \@sff_files,
											version => $version,
											version_subdirectory => $subdir,
											);
    isa_ok($version_run_assembly,'Genome::Model::Tools::454::Newbler::RunAssembly');
    is ($version_run_assembly->newbler_bin, $expected_path, 'found expected path');
}
