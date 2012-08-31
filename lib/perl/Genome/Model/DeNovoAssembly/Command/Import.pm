package Genome::Model::DeNovoAssembly::Command::Import;

use strict;
use warnings;

use Genome;

class Genome::Model::DeNovoAssembly::Command::Import {
    is => 'Command::V2',
    has_input => [
        build => { is => 'Genome::Model::Build::DeNovoAssembly',
            is_output => 1},
    ],
};


sub execute {
    my $self = shift;

    my $build = $self->build;

    my $target_dir = $build->data_directory.'/assembly_files';
    Genome::Sys->create_directory($target_dir);

    $self->status_message('Copying files from '
        . $build->model->import_location . " to $target_dir\n");

    for my $file ($build->required_files, $build->optional_files) {
        my $name = File::Basename::basename($file);
        Genome::Sys->copy_file($file, "$target_dir/$name");
        $self->status_message("Copied $file to $target_dir/$name");
    }
    return 1;
}

1;
