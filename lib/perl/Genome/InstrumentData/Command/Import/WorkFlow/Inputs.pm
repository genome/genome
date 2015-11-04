package Genome::InstrumentData::Command::Import::WorkFlow::Inputs;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Command::Import::WorkFlow::SourceFiles;

class Genome::InstrumentData::Command::Import::WorkFlow::Inputs { 
    is => 'UR::Object',
    id_by => {
        process_id => { is => 'Text', },
        line_number => { is => 'Number', },
    },
    has => {
        analysis_project_id => { is => 'Text', },
        library_id => { is => 'Text', },
        instrument_data_properties => { is => 'HASH', default_value => {}, },
        source_paths => { is => 'ARRAY', },
    },
    has_transient => {
        format => { via => 'source_files', to => 'format', },
        library_name => { via => 'library', to => 'name', },
        sample_name => { via => 'library', to => 'sample_name', },
    },
};

sub create {
    my ($class, %params) = @_;

    my $self = $class->SUPER::create(%params);
    return if not $self;

    for my $requried (qw/ analysis_project_id library_id source_paths /) {
        die $self->error_message("No $requried given to work flow inputs!") if not $self->$requried;
    }

    if ( not $self->instrument_data_properties->{original_data_path} ) {
        $self->instrument_data_properties->{original_data_path} = join(',', $self->source_files->paths);
    }

    if ( $self->process_id ) {
        $self->{instrument_data_properties}->{process_id} = $self->process_id;
    }

    return $self;
}

sub analysis_project {
    return Genome::Config::AnalysisProject->get(id => $_[0]->analysis_project_id);
}

sub library {
    return Genome::Library->get(id => $_[0]->library_id);
}

sub process {
    return Genome::InstrumentData::Command::Import::Process->get(id => $_[0]->process_id);
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
    return Genome::InstrumentData::Command::Import::WorkFlow::SourceFiles->create(paths => $self->source_paths);
}

sub as_hashref {
    my $self = shift;

    my %hash = map { $_ => $self->$_ } (qw/
        analysis_project instrument_data_properties library library_name sample_name
        /);
    $hash{downsample_ratio} = $self->instrument_data_properties->{downsample_ratio};
    $hash{source_paths} = [ $self->source_files->paths ];

    return \%hash;
}

1;

