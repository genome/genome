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
        entity_params => {
            is => 'HASH',
            default_value => {
                individual => {}, sample => {}, library => {},
            },
        },
        source_paths => { is => 'Text', is_many => 1, },
    },
    has_optional => {
        analysis_project_id => { is => 'Text', },
        base_working_directory => { is => 'Text' },
    },
    has_transient => {
        format => { via => 'source_files', to => 'format', },
    },
};

sub get { shift; Genome::InstrumentData::Command::Import::Inputs::Factory->create->from_inputs_id(@_); }
sub create { Carp::confess('Use inputs factory to create!'); }

sub lib_and_source_file_md5sum {
    my $self = shift;
    return substr(
        Genome::Sys->md5sum_data( join(' ', $self->library->name, $self->source_paths) ), 
        0, 6,
    );
}

sub analysis_project {
    my $self = shift;
    if ( $self->process ) {
        return $self->process->analysis_project;
    }
    elsif ( $self->analysis_project_id ) {
        return Genome::Config::AnalysisProject->get(id => $self->analysis_project_id);
    }
    $self->fatal_message('No process or analysis_project_id to get analysis project!');
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
    return Genome::InstrumentData::Command::Import::Inputs::SourceFiles->create(paths => [$self->source_paths]);
}

sub as_hashref {
    my $self = shift;

    my %hash = map { $_ => $self->$_ } (qw/ analysis_project library /);
    $hash{instrument_data_properties} = $self->entity_params->{instdata};
    $hash{downsample_ratio} = $self->entity_params->{instdata}->{downsample_ratio} if $self->entity_params->{instdata}->{downsample_ratio};
    my @source_paths = $self->source_paths;
    $hash{source_paths} = \@source_paths;

    return \%hash;
}

1;

