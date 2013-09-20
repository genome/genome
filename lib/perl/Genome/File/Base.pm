package Genome::File::Base;
use strict;
use warnings;
use Genome;

class Genome::File::Base {
    is => 'UR::Value',
    id_by => [
        path => { is => 'FilesystemPath' }, 
    ],
};

sub open {
    my $self = shift;
    my $type = shift;
    my $path = $self->path;
    # TODO: switch Genome::Sys to use this, not the other way around
    # and have this take any paarm usable in IO::File
    if (not defined $type or $type eq 'r') {
        return Genome::Sys->open_file_for_reading($path);
    }
    elsif($type eq 'w') {
        return Genome::Sys->open_file_for_writing($path);
    }
    else {
        die "unsuppoerted open type $type\n"; 
    }
}

1;

