package Genome::InstrumentData::Command::Import::Process;

use strict;
use warnings FATAL => 'all';

use Genome;

require Cwd;
require File::Basename;
require File::Spec;
require Genome::Sys;

class Genome::InstrumentData::Command::Import::Process {
    is => 'Genome::Process',
    has_input => {
        import_file => { is => 'Text', }, # original location
        import_md5 => { is => 'Text', },
    },
};

sub create {
    my ($class, %params) = @_;

    my $import_file = delete $params{import_file};
    die $class->error_message('No import file given to create process!') if not defined $import_file;
    $import_file = Cwd::abs_path($import_file);
    die $class->error_message('import file (%s) given to create process does not exist!', $import_file) if not -s $import_file;

    my $md5 = Genome::Sys->md5sum($import_file);
    die $class->error_message('No md5 for import file! %s', $import_file) if not $md5;

    return $class->SUPER::create(
        import_file => $import_file,
        import_md5 => $md5,
    );
}

sub create_disk_allocation {
    my $self = shift;

    my $kb_requested = ( -s $self->import_file ) / 1024;
    my $rv = $self->SUPER::create_disk_allocation($kb_requested);
    
    # Copy the import file
    Genome::Sys->copy_file($self->import_file, $self->saved_import_file);
    $self->debug_message("Copied import file to: %s", $rv->absolute_path);

    return $rv;
}

sub saved_import_file { # import_directory location
    my $self = shift;
    return if not $self->disk_allocation;
    my $import_file_basename = File::Basename::basename($self->import_file);
    return File::Spec->join($self->metadata_directory, $import_file_basename);
}

sub instrument_data {
    return map { $_->instrument_data } Genome::InstrumentDataAttribute->get(
        attribute_label => 'process_id',
        attribute_value => $_[0]->id
    );
}

1;

