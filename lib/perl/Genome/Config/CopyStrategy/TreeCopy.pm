package Genome::Config::CopyStrategy::TreeCopy;

use warnings;
use strict;

use Genome;
use File::Spec;
use File::Path 'make_path';
use File::Basename;
use File::Copy;

class Genome::Config::CopyStrategy::TreeCopy {
    is => 'UR::Object',
    doc => 'This copy strategy will mirror a tree of directories with a config file at the bottom between two different roots'
};

sub copy_config {
    my ($self, $file, $source_root, $destination_root) = @_;
    die('You must specify a file, a source, and a destination root!')
        unless $file && $source_root && $destination_root;
    die("$source_root or $destination_root is not a directory!")
        unless (-d $source_root && -d $destination_root);
    die("$file isn't a file!") unless (-f $file);

    _make_directory_structure($file, $source_root, $destination_root);
    _copy_file($file, $source_root, $destination_root);
}

sub _make_directory_structure {
    my ($source_file, $source_root, $dest_root) = @_;
    my $directories_to_create = dirname(_switch_file_root($source_file,
            $source_root, $dest_root));
    make_path($directories_to_create);
}

sub _copy_file {
    my ($source_file, $source_root, $dest_root) = @_;
    copy($source_file, _switch_file_root($source_file, $source_root, $dest_root))
        or die("Failed to copy $source_file to $dest_root!");
}

sub _switch_file_root {
    my ($source_file, $source_root, $dest_root) = @_;
    return $dest_root . '/' .
        File::Spec->abs2rel($source_file, $source_root)
}

1;
