package Genome::Db::Ensembl::VepPlugin;

use strict;
use warnings FATAL => 'all';
use Genome;
use IPC::Run qw(run);

my $INSTALLATION_DIR = '/usr/lib/';

class Genome::Db::Ensembl::VepPlugin {
    has => [
        descriptor => {
            is => 'String',
            doc => 'A string like "Condel@PLUGIN_DIR@b@2"',
        },
        version => {
            is => 'String',
            default => "1",
            doc => 'Version of the vepplugins package to use',
        },
        staging_directory => {
            is => 'Path',
            doc => 'The location where this plugin will stage itself into',
        },
    ],
};

sub stage {
    my $self = shift;

    $self->copy_to_staging_directory;
    $self->replace_configuration_paths;
}

sub copy_to_staging_directory {
    my $self = shift;

    Genome::Sys->create_symlink($self->source_file, $self->dest_file);
    Genome::Sys->rsync_directory(source_directory => $self->source_dir,
        target_directory => $self->dest_dir);
}

sub replace_configuration_paths {
    my $self = shift;

    my $config_dir = File::Spec->join($self->dest_dir, 'config');
    my $pattern = "s|path/to/config|$config_dir|";
    for my $file (glob File::Spec->join($config_dir, '*.conf')) {
        my @sed_cmd = (
            'sed', '-i', $pattern, $file,
        );
        run(@sed_cmd);
    };
}

sub command_line_args {
    my $self = shift;

    return ('--plugin', $self->args);
}

sub args {
    my $self = shift;

    my $config_dir = File::Spec->join($self->dest_dir, 'config');
    (my $processed_args = $self->descriptor) =~ s/PLUGIN_DIR/$config_dir/;
    $processed_args =~ s/\@/,/g;
    return $processed_args;
}

sub name {
    my $self = shift;

    return ($self->_fields)[0];
}

sub _fields {
    my $self = shift;

    return split(/@/, $self->descriptor);
}

sub source_dir {
    my $self = shift;

    return File::Spec->join($self->_installation_dir, 'config', $self->name);
}

sub source_file {
    my $self = shift;

    my $filename = $self->name . '.pm';
    return File::Spec->join($self->_installation_dir, $filename);
}

sub _installation_dir {
    my $self = shift;

    my $subdir = 'vepplugins-' . $self->version;
    return File::Spec->join($INSTALLATION_DIR, $subdir);
}

sub dest_dir {
    my $self = shift;

    return File::Spec->join($self->staging_directory, 'config', $self->name);
}

sub dest_file {
    my $self = shift;

    my $filename = $self->name . '.pm';
    return File::Spec->join($self->staging_directory, $filename);
}

