package Genome::InstrumentData::Command::Import::WorkFlow::SourceFiles;

use strict;
use warnings;

use Genome;

require File::Basename;
require File::Copy;
require Filesys::Df;
use Genome::InstrumentData::Command::Import::WorkFlow::SourceFile;
require List::MoreUtils;
use Params::Validate ':types';

class Genome::InstrumentData::Command::Import::WorkFlow::SourceFiles { 
    is => 'UR::Object',
    has => {
        paths => { is => 'Text', is_many => 1, },
    },
    has_transient => {
        source_files => { is => 'Genome::InstrumentData::Command::Import::WorkFlow::SourceFile', is_many => 1 },
        format => { is =>'Text', },
        retrieval_method => { is =>'Text', },
    },
};

sub original_data_path {
    return join(',', $_[0]->paths);
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    die $self->error_message('No source files given!') if not $self->paths;

    $self->source_files([
        map { 
            Genome::InstrumentData::Command::Import::WorkFlow::SourceFile->create(path => $_)
        } $self->paths
        ]);

    $self->_resolve_property_for_source_files('format');
    $self->_resolve_property_for_source_files('retrieval_method');

    return $self;
}

sub _resolve_property_for_source_files {
    my ($self, $property) = Params::Validate::validate_pos(@_, {type => OBJECT}, {type => SCALAR});

    my @values;
    for my $source_file ( $self->source_files ) {
        push @values, $source_file->$property;
    }

    @values = List::MoreUtils::uniq(@values);
    if ( @values > 1 ) {
        die $self->error_message('Mixed values for source file %s! %s', $property, join(' ', $self->paths));
    }

    return $self->$property($values[0]);
}

sub kilobytes_required_for_processing {
    my $self = shift;

    my $kb_required = 0;
    for my $source_file ( $self->source_files ) {
        $kb_required += $source_file->kilobytes_required_for_processing;
    }

    die $self->error_message('Failed to calculate kilobytes required for processing!') if not $kb_required;

    return $kb_required;
}

sub verify_adequate_disk_space_is_available_for_processing {
    my ($self, $directory) = Params::Validate::validate_pos(@_, {type => OBJECT}, {type => SCALAR});
    $self->debug_message('Verify adequate disk space is available...');

    $self->debug_message("Directory: $directory");
    my $df = eval{ Filesys::Df::df($directory); };
    if( not $df ) {
        $self->error_message($@) if $@;
        die $self->error_message('Failed to get "df" for temp dir! %s', $directory);
    }
    $self->debug_message("Available Kb: ".$df->{bavail});

    my $kb_required = $self->kilobytes_required_for_processing;
    $self->debug_message("Required Kb: ".$kb_required);

    my $remaining_space = $df->{bavail} - $kb_required;
    if ( $remaining_space < 1024 ) { # 1 Mb
        die $self->error_message("There is not enough space in %s to process source files!", $directory);
    }

    $self->debug_message('Verify adequate disk space is available...done');
    return 1;
}

1;

