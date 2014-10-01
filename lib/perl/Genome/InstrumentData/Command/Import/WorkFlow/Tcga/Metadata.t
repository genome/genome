#! /gsc/bin/perl

BEGIN {
    $ENV{UR_COMMAND_DUMP_DEBUG_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Test::Exception;
use Test::More;

my $class = 'Genome::InstrumentData::Command::Import::WorkFlow::Tcga::Metadata';
use_ok($class) or die;

# Failures 
throws_ok(sub{ $class->create(); } , qr/Need UUID or existing metadata file!/, 'create w/o uuid and existing metadata_file');
throws_ok(sub{ $class->create(metadata_file => 'blah'); } , qr/Need UUID or existing metadata file!/, 'create w/o uuid and existing metadata_file');
throws_ok(sub{ $class->create(uuid => 'INVALID'); }, qr/\QFailed to find uuid (INVALID) on CG Hub!\E/, 'create w/ invalid uuid fails');

# Success - uuid
my $uuid = '387c3f70-46e9-4669-80e3-694d450f2919';
my $metadata1 = $class->create(
    uuid => $uuid,
);
ok($metadata1, 'create w/ uuid') or die;
is($metadata1->uuid, $uuid, 'uuid');
ok(-s $metadata1->metadata_file, 'metadata_file downloaded'); # should be temp file
ok($metadata1->_metadata, '_metadata');
ok(_test_metadata($metadata1));

# Success - load
my $metadata2 = $class->create(
    metadata_file => $metadata1->metadata_file,
);
ok($metadata2, 'create w/ metadata_file') or die;
is($metadata2->uuid, $uuid, 'uuid');
is($metadata2->metadata_file, $metadata1->metadata_file, 'metadata_file'); # should be temp file
ok($metadata2->_metadata, '_metadata');
ok(_test_metadata($metadata2));

# Fails to get
throws_ok(sub{ $metadata2->checksum_content_for_file_name(); }, qr/No file name given to get attribute value!/, 'checksum_content_for_file_name fails w/o file name');
throws_ok(sub{ $metadata2->checksum_type_for_file_name(); }, qr/No file name given to get attribute value!/, 'checksum_type_for_file_name fails w/o file name');

done_testing();

###

sub _test_metadata {
    my $metadata = shift;

    # attrs values
    throws_ok(sub{ $metadata->get_attribute_value(); }, qr/No name to get attribute value!/, 'get_attribute_value fails w/o attribute name');
    my %expected_instdata_attrs = (
        aliquot_id => 'f7de2e89-ee90-4098-b86e-57a489b3a71a',
        analysis_id => $uuid,
        center_name => 'BI',
        import_source_name => 'BI',
        library_strategy => 'WXS',
        tcga_name => 'TCGA-77-8154-10A-01D-2244-08',
        sample_id => 'f39b4cc9-9253-4cf9-8827-ebf26af1003a',
        target_region => 'agilent_sureselect_exome_version_2_broad_refseq_cds_only_hs37',
        uuid => $uuid,
    );
    for my $name ( sort keys %expected_instdata_attrs ) {
        is($metadata->get_attribute_value($name), $expected_instdata_attrs{$name}, "get_attribute_value for $name");
    }

    # files
    my @expected_files = (
        { 
            file_name => 'C508.TCGA-77-8154-10A-01D-2244-08.1.bam',
            file_size => 12322789137,
            checksum_content => '5d2c5cbfc7420405fd4e8e7491a56dc8',
        },
        {
            file_name => 'C508.TCGA-77-8154-10A-01D-2244-08.1.bam.bai',
            file_size => 5731680,
            checksum_content => '10b2dd0951d76c4813e91e33074d1680',
        },
    );
    is_deeply([$metadata->file_names], [map { $_->{file_name} } @expected_files], 'file_names');
    is_deeply([$metadata->bam_file_names], [$expected_files[0]->{file_name}], 'bam_file_names');

    # file attrs
    for my $file ( @expected_files ) {
        my $file_name = $file->{file_name};
        is($metadata->checksum_type_for_file_name($file_name), 'MD5', "checksum_content_for_file_name for $file_name");
        is(
            $metadata->checksum_content_for_file_name($file_name),
            $file->{checksum_content},
            "checksum_type_for_file_name for $file_name",
        );
        is(
            $metadata->filesize_in_kb_for_file_name($file_name),
            $file->{file_size},
            "checksum_type_for_file_name for $file_name",
        );
    }

    return 1;
}

