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

    $self->debug_message('Copying files from '
        . $build->model->import_location . " to $target_dir\n");

    for my $name ( $self->_required_file_names ) {
        my $file = $build->model->import_location."/$name";
        if ( not -e $file ) {
            $self->warning_message("Required file name, $name, is missing: $file");
            next;
        }
        Genome::Sys->copy_file( $file, "$target_dir/$name");
        $self->debug_message("Copied $file to $target_dir");
    }

    for my $name ( $self->_optional_file_names ) {
        my $file = $build->model->import_location."/$name";
        if ( -e $file ) {
            Genome::Sys->copy_file( $file, "$target_dir/$name");
            $self->debug_message("Copied $file to $target_dir");
        }
    }

    return 1;
}


sub _required_file_names {
    return qw/
ASSEMBLER
AUTHOR
BIOPROJECT
COVERAGE
READ_TYPE
RELEASE_NOTES
contigs.bases
supercontigs.fasta
supercontigs.agp
/;
}

sub _optional_file_names {
    return qw/
BUILD_ID
README
readme
VERSION
version
stats.txt
basic.stats.txt
contigs.quals
supercontigs.quals
/;
}

1;
