package Genome::Model::Build::DeNovoAssembly::Import;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::DeNovoAssembly::Import {
    is => 'Genome::Model::Build::DeNovoAssembly',
};

sub required_files {
    my $self = shift;

    my @required_files;
    for my $name ( $self->_required_file_names ) {
        my $file = $self->model->import_location."/$name";
        if ( not -s $file ) {
            die "Required file, $name, is missing for empty: $file";
        }
        push @required_files, $file;
    }
    return @required_files;
}

sub optional_files {
    my $self = shift;

    my @optional_files;
    for my $name ( $self->optional_file_names ) {
        my $file = $self->model->import_location."/$name";
        push @optional_files, $file if -s $file;
    }
    return @optional_files;
}

sub _required_file_names {
    return qw/
contigs.bases
supercontigs.fasta
supercontigs.agp
/;
}

sub optional_file_names {
    return qw/
README
readme
contigs.quals
supercontigs.quals
/;
}

1;
