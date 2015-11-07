package Genome::InstrumentData::Command::Import::Inputs;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Command::Import::Inputs::SourceFiles;

class Genome::InstrumentData::Command::Import::Inputs {
    is => 'UR::Object',
    id_by => [ # use arrayref to keep order
        process_id => { is => 'Text', },
        line_number => { is => 'Number', },
    ],
    has => {
        analysis_project_id => { is => 'Text', },
        entity_params => {
            is => 'HASH',
            default_value => {
                individual => {}, sample => {}, library => {},
            },
        },
        source_paths => { is => 'ARRAY', },
    },
    has_transient => {
        format => { via => 'source_files', to => 'format', },
        library_name => { via => 'library', to => 'name', },
        sample_name => { via => 'library', to => 'sample_name', },
    },
};

sub lib_and_source_file_md5sum {
    my $self = shift;
    return substr(
        Genome::Sys->md5sum_data( join(' ', $self->library_name, @{$self->source_paths}) ), 
        0, 6,
    );
}

sub create {
    my ($class, %params) = @_;

    my $self = $class->SUPER::create(%params);
    return if not $self;

    for my $requried (qw/ analysis_project_id source_paths /) {
        die $self->error_message("No $requried given to work flow inputs!") if not $self->$requried;
    }

    if ( not $self->entity_params->{instdata}->{original_data_path} ) {
        $self->entity_params->{instdata}->{original_data_path} = join(',', $self->source_files->paths);
    }

    if ( $self->process_id ) {
        $self->entity_params->{instdata}->{process_id} = $self->process_id;
    }

    return $self;
}

sub analysis_project {
    return Genome::Config::AnalysisProject->get(id => $_[0]->analysis_project_id);
}

sub library {
    my $self = shift;
    my %params;
    if ( $self->entity_params->{library}->{id} ) {
        %params = ( id => $self->entity_params->{library}->{id} );
    }
    elsif ( $self->entity_params->{library}->{name} ) {
        %params = ( name => $self->entity_params->{library}->{name} );
    }
    else {
        $self->fatal_message('No library id or name given to inputs!');
    }
    return Genome::Library->get(%params);
}

sub process {
    return Genome::InstrumentData::Command::Import::Process->get(id => $_[0]->process_id);
}

sub instrument_data_properties {
    return $_[0]->entity_params->{instdata};
}

sub instrument_data_for_original_data_path {
    my $self = shift;
    my @odp_attrs = Genome::InstrumentDataAttribute->get(
        attribute_label => 'original_data_path',
        attribute_value => $self->source_files->original_data_path,
    );
    return if not @odp_attrs;
    return map { $_->instrument_data } @odp_attrs;
}

sub source_files {
    my $self = shift;
    return Genome::InstrumentData::Command::Import::Inputs::SourceFiles->create(paths => $self->source_paths);
}

sub as_hashref {
    my $self = shift;

    my %hash = map { $_ => $self->$_ } (qw/
        analysis_project library library_name sample_name
        /);
    $hash{instrument_data_properties} = $self->entity_params->{instdata};
    $hash{downsample_ratio} = $self->entity_params->{instdata}->{downsample_ratio};
    $hash{source_paths} = [ $self->source_files->paths ];

    return \%hash;
}

1;

