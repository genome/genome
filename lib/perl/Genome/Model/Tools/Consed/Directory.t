#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 5;
use Test::MockObject;

use Genome::Model::Tools::Consed::Directory;
use File::Path;

=cut
use Test::Mock;
# Mock up a processing profile
my $processing_profile = Test::MockObject->new();
$processing_profile->fake_module('Genome::ProcessingProfile');
$processing_profile->set_always('type_name', 'reference alignment');
=cut

my $path = Genome::Sys->create_temp_directory;
rmtree $path; # Want the path, but not the directory... yet

my $consed_dir = Genome::Model::Tools::Consed::Directory->create(directory => $path);
ok (!$consed_dir, "Didn't create a Genome::Model::Tools::Consed::Directory without an existing directory");

system "touch $path";
$consed_dir = Genome::Model::Tools::Consed::Directory->create(directory => $path);
ok (!$consed_dir, "Didn't create a Genome::Model::Tools::Consed::Directory with a non-directory");

unlink $path;

create_test_fixture();

$consed_dir = Genome::Model::Tools::Consed::Directory->create(directory => $path);
ok ($consed_dir, "created a Genome::Model::Tools::Consed::Directory");


my @directories = $consed_dir->directories;
ok ($directories[0] eq 'edit_dir' && $directories[1] eq 'phd_dir' && $directories[2] eq 'chromat_dir', "Got correct directories");

$consed_dir->create_consed_directory_structure;
ok(-d "$path/".$directories[0] && -d "$path/".$directories[0] && -d "$path/".$directories[0], "Created directory structure");

destroy_test_fixture();

sub create_test_fixture {
    mkdir $path;
}

sub destroy_test_fixture {
    rmtree $path;    
}
