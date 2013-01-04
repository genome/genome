#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
if (Genome::Config->arch_os ne 'x86_64') {
   plan skip_all => 'requires 64-bit machine';
}

use above 'Genome';
use Data::Dumper 'Dumper';

use_ok('Genome::InstrumentData::Solexa::Report::Quality') or die;

my $base_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-InstrumentData-Solexa-Report-Quality/';
ok( -d $base_dir, "Base dir exists" ) or die;

my $bam = $base_dir.'/test_run_name.sanger.bam';
ok( -s $bam, "Bam file exists" ) or die;

my $archive = $base_dir.'/test_run_name.sanger.tar.gz';
ok( -s $archive, "Archive file exists" ) or die;

my $ori_report_xml   = $base_dir .'/report.xml';
ok( -s $ori_report_xml, "Report xml exists" );

my $tmp = Genome::Sys->create_temp_directory();

my @test_types = ( qw/ bam archive /);
my %input = (
    bam_path => $bam,
    archive_path => $archive,
);

for my $input_type ( qw/ bam_path archive_path / ) {

    my $sub_dir_name = shift @test_types;

    my $sub_dir = Genome::Sys->create_directory( $tmp."/$sub_dir_name" );
    my %params = &inst_data_params;
    $params{$input_type} = $input{$input_type};

    my $inst_data = Genome::InstrumentData::Solexa->create_mock( %params );
    isa_ok( $inst_data, 'Genome::InstrumentData::Solexa' );

    my $r = Genome::InstrumentData::Solexa::Report::Quality->create(
        instrument_data_id => $inst_data->id,
    );
    ok($r, "created a new report");
    
    my $v = $r->generate_report;
    ok($v, "generation worked");

    my $result = $v->save($sub_dir);
    ok($result, "saved to $sub_dir");
    
    my $name = $r->name;
    $name =~ s/ /_/g;
    
    ok(-d "$sub_dir/$name", "report directory $sub_dir/$name is present");
    ok(-e "$sub_dir/$name/report.xml", 'xml report is present');
    
    my @diff = `diff "$sub_dir/$name/report.xml" $ori_report_xml`;
    is(scalar @diff, 4, 'report.xml is created as expected'); #Only time stamp is different
}

#<STDIN>;

done_testing();



sub inst_data_params {
    return (
        id => '-123456',
        sequencing_platform => 'solexa',
        sample_name => 'test_sample_name',
        library_name => 'test_library_name',
        library_id => '-1233445',
        run_name => 'test_run_name',
        subset_name => 4,
        lane => 4,
        run_type => 'Paired End Read 2',
        is_paired_end => 1,
    );
}
