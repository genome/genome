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
        analysis_project => { is => 'Genome::Config::AnalysisProject', },
        import_file => { is => 'Text', }, # original location
        import_md5 => { is => 'Text', },
    },
};

sub create {
    my ($class, %params) = @_;

    die $class->error_message('No import file given to create process!') if not $params{import_file};
    $params{import_file} = Cwd::abs_path($params{import_file});
    die $class->error_message('Import file (%s) given to create process does not exist!', $params{import_file}) if not -s $params{import_file};

    $params{import_md5} = Genome::Sys->md5sum($params{import_file});
    die $class->error_message('No md5 for import file! %s', $params{import_file}) if not $params{import_md5};

    return $class->SUPER::create(%params);
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

