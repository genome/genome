#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
};

use above "Genome";

require Genome::Utility::Test;
require File::Compare;
require File::Spec;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::Basic') or die;

# Test data
my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', File::Spec->catfile('cghub', 'v2')) or die;
my $source_bam = File::Spec->catfile($data_dir, '387c3f70-46e9-4669-80e3-694d450f2919.bam');

# Create library
use_ok('Genome::Sample::Command::Import') or die;
my $sample_importer = Genome::Sample::Command::Import::Tcga->execute(
    name => 'TCGA-AB-2804-03B-01W-0728-08',
);
my $library = $sample_importer->_library;
ok($library, 'created library') or die;

# Local Bam
my $cmd = Genome::InstrumentData::Command::Import::Basic->execute(
    library => $library,
    import_source_name  => 'CGHub',
    source_files => [$source_bam],
    instrument_data_properties => [qw/ target_region_set_name=none tcga_name=TCGA-AB-2804-03B-01W-0728-08 /],
);
ok($cmd->result, "execute");

my @instdata =  $cmd->_new_instrument_data;
is(@instdata, 1, 'create instrument data');
my $instdata = $instdata[0];

# Bam
my $bam_path = $instdata->bam_path;
ok(-s $bam_path, "bam exists");
my $expected_bam_path = File::Spec->join($data_dir, 'expected.bam');
is(File::Compare::compare($bam_path, $expected_bam_path), 0, 'bam matches');
is(File::Compare::compare($bam_path.'.flagstat', $expected_bam_path.'.flagstat'), 0, 'flagstat matches');

# Attrs
is($instdata->import_source_name, 'BI', 'import source name'); # retrieved from metadata
is($instdata->sequencing_platform, 'solexa', 'platform');
is($instdata->original_data_path, $source_bam, 'original data path');
is($instdata->user_name, Genome::Sys->username, "user name is correct");
is($instdata->target_region_set_name, 'none', "target_region_set_name is correct");
ok($instdata->import_date, "date is set");

# TCGA attrs
my %expected_attributes = (
    aliquot_id => 'f7de2e89-ee90-4098-b86e-57a489b3a71a',
    analysis_id => '387c3f70-46e9-4669-80e3-694d450f2919',
    participant_id => '569691a3-15b4-4b1c-b8b7-b3ad17d0996e',
    sample_id => 'f39b4cc9-9253-4cf9-8827-ebf26af1003a',
    tcga_name => 'TCGA-AB-2804-03B-01W-0728-08',
);
for my $label ( keys %expected_attributes ) {
    is(
        eval{ $instdata->attributes(attribute_label => $label)->attribute_value },
        $expected_attributes{$label},
        "$label has correct value"
    );
}

# Metrics
my ($expected_read_count, $expected_read_length) = (qw/ 996 35 /);
is($instdata->read_count, $expected_read_count, "read_count is set");
is($instdata->read_length, $expected_read_length, "read_length is correct");
ok($instdata->is_paired_end, "is_paired_end is correct");

done_testing();
