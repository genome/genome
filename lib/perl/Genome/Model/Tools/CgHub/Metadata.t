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

my $metadata = $class->create;
ok($metadata, 'create metadata');

# Failures 
throws_ok(sub{ $metadata->add_xml; }, qr/add_xml but 2 were expected/, 'create failed w/o XML');
throws_ok(sub{ $metadata->add_xml('blah'); }, qr//, 'create failed w/ invalid xml');
throws_ok(sub{ $metadata->add_file('blah'); }, qr/File \(blah\) does not exist/, 'add_file failed w/ non existing file');

# Success [b37]
my $analysis_id1 = '387c3f70-46e9-4669-80e3-694d450f2919';
ok($metadata, 'create') or die;
ok($metadata->add_file($xml_file), 'add_file to metadata') or die;
ok($metadata->metadata, 'metadata');

# Fails to get_attribute_value
throws_ok(sub{ $metadata->get_attribute_value(); }, qr/but 3 were expected/, 'get_attribute_value fails w/o attribute name');
throws_ok(sub{ $metadata->get_attribute_value($analysis_id1); }, qr/but 3 were expected/, 'get_attribute_value fails w/o attribute name');
# Fails to get checksum type/content
throws_ok(sub{ $metadata->checksum_content_for_file_name(); }, qr/but 3 were expected/, 'checksum_content_for_file_name fails w/o file name');
throws_ok(sub{ $metadata->checksum_type_for_file_name($analysis_id1); }, qr/but 3 were expected/, 'checksum_type_for_file_name fails w/o file name');

ok(
    _test_metadata(
        metadata => $metadata,
        lookup => {
            legacy_sample_id => 'TCGA-77-8154-10A-01D-2244-08',
            analysis_id => $analysis_id1,
        },
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
            center_name => 'BI',
            library_strategy => 'WXS',
            sample_id => 'f39b4cc9-9253-4cf9-8827-ebf26af1003a',
        },
    ),
    'test metadata'
);

done_testing();

###

sub _test_metadata {
    my %params = @_;

    my $metadata = $params{metadata};
    my $analysis_id = $params{lookup}->{analysis_id};
    my $legacy_sample_id = $params{lookup}->{legacy_sample_id};

    for my $lookup1 (qw/ analysis_id legacy_sample_id /) {
        for my $lookup2 (qw/ analysis_id legacy_sample_id /) {
            is(
                $metadata->get_attribute_value($params{lookup}->{$lookup1}, $lookup1),
                $params{lookup}->{$lookup1},
                "get_attribute_value for $params{lookup}->{$lookup1} $lookup1",
            );
            is(
                $metadata->get_attribute_value($params{lookup}->{$lookup1}, $lookup2),
                $params{lookup}->{$lookup2},
                "get_attribute_value for $params{lookup}->{$lookup1} $lookup2",
            );
        }
    }

    # attrs values
    my %instdata_attrs = %{$params{instdata_attrs}};
    for my $name ( sort keys %instdata_attrs ) {
        is($metadata->get_attribute_value($analysis_id, $name), $instdata_attrs{$name}, "get_attribute_value for $analysis_id $name");
        is($metadata->get_attribute_value($legacy_sample_id, $name), $instdata_attrs{$name}, "get_attribute_value for $legacy_sample_id $name");
    }

    # files
    my @files = @{$params{files}};
    my @file_names = $metadata->file_names($analysis_id);
    is_deeply(\@file_names, [map { $_->{file_name} } @files], 'file_names');
    is_deeply([$metadata->bam_file_names($analysis_id)], [$files[0]->{file_name}], 'bam_file_names');

    # file attrs
    for my $file_name ( @file_names ) {
        my ($file) = grep { $_->{file_name} eq $file_name } @files;
        ok($file, "found file hash for $file_name");
        is(
            $metadata->checksum_type_for_file_name($analysis_id, $file_name),
            $file->{checksum}->{type},
            "checksum_content_for_file_name for $analysis_id $file_name",
        );
        is(
            $metadata->checksum_content_for_file_name($analysis_id, $file_name),
            $file->{checksum}->{content},
            "checksum_type_for_file_name for $analysis_id $file_name",
        );
        is(
            $metadata->file_size_for_file_name($analysis_id, $file_name),
            $file->{file_size},
            "checksum_type_for_file_name for $analysis_id $file_name",
        );
    }

    return 1;
}

