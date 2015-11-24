package Genome::SoftwareResult::StageableSimple::SingleFile;

use strict;
use warnings;
use Genome;
use File::Basename;

class Genome::SoftwareResult::StageableSimple::SingleFile {
    is => 'Genome::SoftwareResult::StageableSimple',
    is_abstract => 1,
};

sub _temp_staging_file_path {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, $self->_file_name);
}

sub _file_name {
    my $self = shift;
    Carp::confess("Abstract method (_file_name). Needs to be overwritten in your child class");
}

sub file_path {
    my $self = shift;
    return File::Spec->join($self->output_dir, $self->_file_name);
}

1;
