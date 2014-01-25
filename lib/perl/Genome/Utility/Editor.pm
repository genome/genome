package Genome::Utility::Editor;

use strict;
use warnings;

use Genome;

sub from_empty_file {
    my $file_path = Genome::Sys->create_temp_file_path();
    return _invoke_editor($file_path, @_);
}

sub from_existing_file {
    my $file_path = shift;
    die("Must give a file name!") unless $file_path;
    my $contents = Genome::Sys->read_file($file_path);

    return from_existing_contents($contents, @_);
}

sub from_existing_contents {
    my $contents = shift;
    die("Must supply file contents!") unless $contents;

    my $file_path = Genome::Sys->create_temp_file_path();
    Genome::Sys->write_file($file_path, $contents);

    return _invoke_editor($file_path, @_);
}

sub _invoke_editor {
    my $file = shift;
    my %params = @_;

    my $allow_empty = delete $params{allow_empty} || 0;

    system(_editor(), $file);

    return undef unless (-e $file);
    return undef unless (-s $file || $allow_empty);
    return Genome::Sys->read_file($file);
}

sub _editor {
    return $ENV{EDITOR} || 'vim';
}

1;
