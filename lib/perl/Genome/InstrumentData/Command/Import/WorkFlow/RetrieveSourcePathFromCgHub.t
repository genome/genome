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
my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', File::Spec->catfile('cghub', 'v1')) or die;
my $working_directory = File::Temp::tempdir(CLEANUP => 1);

my $uuid = '387c3f70-46e9-4669-80e3-694d450f2919';
my $source_basename = $uuid.'.bam';
my $query = Genome::Model::Tools::CgHub::Query->create(uuid => $uuid);
Sub::Install::reinstall_sub({ # do not execute query
        code => sub {
            Genome::Sys->create_symlink(
                File::Spec->join($data_dir, 'metadata.xml'),
                File::Spec->join($working_directory, 'metadata.xml'),
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
            Genome::Sys->create_symlink(
                File::Spec->join($data_dir, $source_basename),
                File::Spec->join($working_directory, $source_basename),
            );
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

my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromCgHub->execute(
    working_directory => $working_directory,
    source_path => $gene_torrent->source_url,
);
ok($cmd->result, 'execute');
my $destination_path = $cmd->destination_path;
is($destination_path, File::Spec->join($working_directory, $source_basename), 'destination path named correctly');
ok(-s $destination_path, 'destination path exists');

my $destination_md5_path = $cmd->destination_md5_path;
my $original_md5_path = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->original_md5_path_for($destination_path);
is($destination_md5_path, $original_md5_path, 'correcly named destination md5 path');
ok(-s $destination_md5_path, 'destination md5 path exists');
my $md5 = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->load_md5($destination_md5_path);
is($md5, 'f81fbc3d3a6b57d11e60b016bb2c950c', 'correct md5');

#print "$working_directory\n"; <STDIN>;
done_testing();
