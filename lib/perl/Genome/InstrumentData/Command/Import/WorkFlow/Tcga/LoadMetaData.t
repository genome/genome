#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::Tcga::LoadMetaData') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'bam-tcga/v1') or die;
my $metadata_path = $test_dir.'/metadata.xml';

# Failures
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::Tcga::LoadMetaData->create(
    metadata_path => 'blah',
);
my @errors = eval{ $cmd->__errors__; };
is(@errors, 1, 'correct number of errors for invalid path');
is($errors[0]->desc, 'Metadata path is not a valid file! blah', 'correct error');

# Success
$cmd = Genome::InstrumentData::Command::Import::WorkFlow::Tcga::LoadMetaData->create(
    metadata_path => $metadata_path,
);
@errors = $cmd->__errors__;
is(@errors, 0, 'no errors');
ok($cmd->execute, 'execute');
my %expected_instdata_attrs = (
    'tcga_name' => 'test_metadata_tcga_name',
    'target_region' => 'agilent_sureselect_exome_version_2_broad_refseq_cds_only_hs37',
    'library_strategy' => 'WXS',
    'bam_md5' => '0f307401916947ab16e37b225da8c919',
    'aliquot_id' => 'f957194b-6da9-4690-a87d-0051e239bf3f',
    'sample_id' => '8aca008c-f55a-420a-82c7-acd2cca77d85',
    'analysis_id' => 'a1d11d67-4d5f-4db9-a61d-a0279c3c3d4f',
);
is_deeply(\%expected_instdata_attrs, $cmd->instrument_data_attributes, 'instrument_data_attributes');

done_testing();
