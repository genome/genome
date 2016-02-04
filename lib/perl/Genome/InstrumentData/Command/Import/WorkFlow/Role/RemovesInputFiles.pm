package Genome::InstrumentData::Command::Import::WorkFlow::Role::RemovesInputFiles;

use strict;
use warnings;

require File::Basename;

use Genome;
use UR::Role 'after';

role Genome::InstrumentData::Command::Import::WorkFlow::Role::RemovesInputFiles {
    has => {
       remove_input_files => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            doc => 'If true, for properties marked as inputs and files, remove their associated values.',
        },
    },
    requires => [qw/ working_directory /],
};

after 'execute' => sub{
    my ($rv, $self, @params) = @_;
    return $rv if not $rv; # execute failed...

    for my $input_file_to_remove ( $self->input_files_to_remove ) {
        unlink $input_file_to_remove;
    }

    return $rv;
};

sub input_file_property_names {
    map { $_->property_name } $_[0]->__meta__->properties(is_input => 1, data_type => 'File');
}

sub input_files_to_remove {
    my $self = shift;

    my @input_files_to_remove;
    for my $input_file_property_name ( $self->input_file_property_names ) {
        my $path = $self->$input_file_property_name;
        # Only remove input files that are defined and are files
        next if not $path or not -f $path;
        # Only remove input files that reside in the working directory
        next if File::Basename::dirname($path) ne $self->working_directory;
        push @input_files_to_remove, glob($path.'*');
    }

    return @input_files_to_remove;
}

1;

