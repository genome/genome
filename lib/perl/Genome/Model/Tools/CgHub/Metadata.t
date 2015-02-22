#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_COMMAND_DUMP_DEBUG_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

require File::Spec;
require Genome::Utility::Test;
use Test::Exception;
use Test::More;

my $class = 'Genome::Model::Tools::CgHub::Metadata';
use_ok($class) or die;

my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::Model::Tools::CgHub');
my $xml_file = File::Spec->join($data_dir, 'metadata.xml');
my $b36_xml_file = File::Spec->join($data_dir, 'metadata.b36.xml');

# Failures 
throws_ok(sub{ $class->create(metadata_file => 'blah'); }, qr/File \(blah\) does not exist/, 'create failed w/ invalid metadata file');

# Success [b37]
my $uuid = '387c3f70-46e9-4669-80e3-694d450f2919';
my $metadata1 = $class->create(
    metadata_file => $xml_file,
);
ok($metadata1, 'create') or die;
ok($metadata1->_metadata, '_metadata');
is($metadata1->uuid, $uuid, 'uuid');
ok(
    _test_metadata(
        metadata => $metadata1,
        files => [
            { 
                file_name => 'C508.TCGA-77-8154-10A-01D-2244-08.1.bam',
                file_size => 12322789137,
                checksum => {
                    content => '5d2c5cbfc7420405fd4e8e7491a56dc8',
                    type => 'MD5',
                },
            },
            {
                file_name => 'C508.TCGA-77-8154-10A-01D-2244-08.1.bam.bai',
                file_size => 5731680,
                checksum => {
                    content => '10b2dd0951d76c4813e91e33074d1680',
                    type => 'MD5',
                },
            },
        ],
        instdata_attrs => {
            aliquot_id => 'f7de2e89-ee90-4098-b86e-57a489b3a71a',
            analysis_id => $uuid,
            center_name => 'BI',
            import_source_name => 'BI',
            library_strategy => 'WXS',
            tcga_name => 'TCGA-77-8154-10A-01D-2244-08',
            sample_id => 'f39b4cc9-9253-4cf9-8827-ebf26af1003a',
            uuid => $uuid,
        },
        target_region => 'agilent_sureselect_exome_version_2_broad_refseq_cds_only_hs37',
        reference_assembly_attrs => {
            shortname => 'HG19_Broad_variant',
            version => '37',
        }
    ),
    'test metadata'
);

# Success - load [b36]
my $uuid2 = 'a1d11d67-4d5f-4db9-a61d-a0279c3c3d4f';
my $metadata2 = $class->create(
    metadata_file => $b36_xml_file,
);
ok($metadata2, 'create w/ b36 xml file') or die;
ok($metadata2->_metadata, '_metadata');
is($metadata2->uuid, $uuid2, 'uuid');
ok(
    _test_metadata(
        metadata => $metadata2,
        files => [
            {
                file_name => 'TCGA-A6-2674-01A-02W-0831-10_SOLiD.bam',
				file_size => 18604187778,
				checksum => {
                    type => "MD5",
                    content => '0f307401916947ab16e37b225da8c919',
                },
            },
            {
                file_name => 'TCGA-A6-2674-01A-02W-0831-10_SOLiD.bam.bai',
                file_size => 7121464,
                checksum => {
                    type => "MD5",
                    content => '89a44bbea9fcd62e5c5e1a3dc3610014',
                },
            },
        ],
        instdata_attrs => {
            aliquot_id => 'f957194b-6da9-4690-a87d-0051e239bf3f',
            analysis_id => $uuid2,
            center_name => 'BCM',
            import_source_name => 'BCM',
            library_strategy => 'WXS',
            tcga_name => 'TCGA-A6-2674-01A-02W-0831-10',
            sample_id => '8aca008c-f55a-420a-82c7-acd2cca77d85',
            uuid => $uuid2,
        },
        target_region => 'agilent sureselect exome version 2 broad refseq cds only',
        reference_assembly_attrs => {
            shortname => 'NCBI36_BCM_variant',
            version => '36',
        }
    ),
    'test metadata'
);

# Fails to get checksum type/content
throws_ok(sub{ $metadata2->checksum_content_for_file_name(); }, qr/No file name given to get attribute value!/, 'checksum_content_for_file_name fails w/o file name');
throws_ok(sub{ $metadata2->checksum_type_for_file_name(); }, qr/No file name given to get attribute value!/, 'checksum_type_for_file_name fails w/o file name');

# Fail to get target region when no library strategy is given
my $library_strategy = $metadata2->_metadata->{Result}->{library_strategy};
$metadata2->_metadata->{Result}->{library_strategy} = undef;
throws_ok(sub{ $metadata2->target_region; }, qr/No library strategy in metadata to resolve target region!/, 'failed to get target region w/o library strategy');
$metadata2->_metadata->{Result}->{library_strategy} = $library_strategy;

done_testing();

###

sub _test_metadata {
    my %params = @_;

    my $metadata = $params{metadata};
    my %instdata_attrs = %{$params{instdata_attrs}};

    # attrs values
    throws_ok(sub{ $metadata->get_attribute_value(); }, qr/No name to get attribute value!/, 'get_attribute_value fails w/o attribute name');
    for my $name ( sort keys %instdata_attrs ) {
        is($metadata->get_attribute_value($name), $instdata_attrs{$name}, "get_attribute_value for $name");
    }

    # ref
    my %reference_assembly_attrs = %{$params{reference_assembly_attrs}};
    is($metadata->reference_assembly_shortname, $reference_assembly_attrs{shortname}, 'reference_assembly_shortname');
    is($metadata->reference_assembly_version, $reference_assembly_attrs{version}, 'reference_assembly_version');

    # files
    my @files = @{$params{files}};
    my @file_names = $metadata->file_names;
    is_deeply(\@file_names, [map { $_->{file_name} } @files], 'file_names');
    is_deeply([$metadata->bam_file_names], [$files[0]->{file_name}], 'bam_file_names');

    # file attrs
    for my $file_name ( @file_names ) {
        my ($file) = grep { $_->{file_name} eq $file_name } @files;
        ok($file, "found file hash for $file_name");
        is(
            $metadata->checksum_type_for_file_name($file_name),
            $file->{checksum}->{type},
            "checksum_content_for_file_name for $file_name",
        );
        is(
            $metadata->checksum_content_for_file_name($file_name),
            $file->{checksum}->{content},
            "checksum_type_for_file_name for $file_name",
        );
        is(
            $metadata->file_size_for_file_name($file_name),
            $file->{file_size},
            "checksum_type_for_file_name for $file_name",
        );
    }

    # target region
    is($metadata->target_region, $params{target_region}, 'correct target_region');

    return 1;
}

