#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Spec;
require File::Temp;
use Test::More;

my $class = 'Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromCgHub';
use_ok($class) or die;
my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', File::Spec->catfile('tcga', 'v3')) or die;
my $working_directory = File::Temp::tempdir(CLEANUP => 1);

use_ok('Genome::Model::Tools::CgHub::GeneTorrent');
use_ok('Genome::Model::Tools::CgHub::Query');
my $uuid = '387c3f70-46e9-4669-80e3-694d450f2919';
my $query = Genome::Model::Tools::CgHub::Query->create(uuid => $uuid);
Sub::Install::reinstall_sub({ # do not execute query
        code => sub {
            my $basename = 'metadata.xml';
            Genome::Sys->create_symlink(
                File::Spec->join($data_dir, $basename),
                File::Spec->join($working_directory, $basename),
            );
            $query->result(1);
            return $query;
        },
        into => $query->class,
        as   => 'execute',
    });

my $gene_torrent = Genome::Model::Tools::CgHub::GeneTorrent->create(uuid => $uuid);
Sub::Install::reinstall_sub({ # do not execute gene torrent
        code => sub {
            my @basenames = (qw/
                387c3f70-46e9-4669-80e3-694d450f2919.bam   
                387c3f70-46e9-4669-80e3-694d450f2919.bam.md5
            /);
            for my $basename ( @basenames ) {
                Genome::Sys->create_symlink(
                    File::Spec->join($data_dir, $basename),
                    File::Spec->join($working_directory, $basename),
                );
            }
            $gene_torrent->result(1);
            return $gene_torrent;
        },
        into => $gene_torrent->class,
        as   => 'execute',
    });
is(
    $class->__meta__->property_meta_for_name('lsf_resource')->default_value,
    $gene_torrent->__meta__->property_meta_for_name('lsf_resource')->default_value,
    'lsf_resource',
);

my $retrieve_from_cghub = Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromCgHub->execute(
    working_directory => $working_directory,
    source_path => $gene_torrent->source_url,
);
ok($retrieve_from_cghub->result, 'execute');
is($retrieve_from_cghub->uuid, $uuid, 'uuid');
is($retrieve_from_cghub->metadata_file, File::Spec->join($working_directory, 'metadata.xml'), 'metadata_file');
my $destination_path = $retrieve_from_cghub->destination_path;
is(
    $destination_path,
    File::Spec->join($working_directory, $uuid.'.bam'),
    'destination path named correctly',
);
ok(-s $destination_path, 'destination path exists');
my $original_md5_path = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->original_md5_path_for($destination_path);
is($original_md5_path, File::Spec->join($working_directory, $uuid.'.bam.md5-orig'), 'destination md5 path named correctly');
ok(-s $original_md5_path, 'destination md5 exists');

#print "$working_directory\n"; <STDIN>;
done_testing();
