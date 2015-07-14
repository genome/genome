package Genome::Model::Build::ReferenceSequence::IndexBase;

use strict;
use warnings;

use Carp;
use Genome;

class Genome::Model::Build::ReferenceSequence::IndexBase {
    is_abstract => 1,
    is => ['Genome::SoftwareResult::Stageable', 'Genome::SoftwareResult::WithNestedResults'],
    has => [
        reference_build         => {
                                    is => 'Genome::Model::Build::ImportedReferenceSequence',
                                    id_by => 'reference_build_id',
                                },
        reference_name          => { via => 'reference_build', to => 'name', is_mutable => 0, is_optional => 1 },
        aligner                 => {
                                    calculate_from => [qw/aligner_name aligner_version aligner_params/],
                                    calculate => q|no warnings; "$aligner_name $aligner_version $aligner_params"|
                                },
    ],
    has_input => [
        reference_build_id      => {
                                    is => 'Number',
                                    doc => 'the reference to use by id',
                                },
    ],
    has_param => [
        aligner_name            => {
                                    is => 'Text', default_value => 'maq',
                                    doc => 'the name of the aligner to use, maq, blat, newbler etc.',
                                },
        aligner_version         => {
                                    is => 'Text',
                                    doc => 'the version of the aligner to use, i.e. 0.6.8, 0.7.1, etc.',
                                    is_optional=>1,
                                },
        aligner_params          => {
                                    is => 'Text',
                                    is_optional=>1,
                                    doc => 'any additional params for the aligner in a single string',
                                },
    ],
    has_transient => [
        aligner_class_name      => {
                                    is => 'Text',
                                    is_optional => 1,
        }
    ]
};

sub full_consensus_path {
    my $self = shift;
    my $extension = shift;

    my $path = $self->output_dir . '/all_sequences';
    if ($extension) {
        $path .= '.'. $extension;
    }
    return $path;
}

sub required_rusage {
    # override in subclasses
    # e.x.: "-R 'span[hosts=1] rusage[tmp=50000:mem=12000]' -M 12000000";
    ''
}

sub resolve_allocation_disk_group_name {
    if ($_[0]->reference_build->model->is_rederivable) {
        return Genome::Config::get('disk_group_models');
    } else {
        return Genome::Config::get('disk_group_references');
    }
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $aligner_name_tag = $self->aligner_name;
    $aligner_name_tag =~ s/[^\w]/_/g;

    my $staged_basename = File::Basename::basename($self->temp_staging_directory);
    my @path_components = ('model_data', $self->_resolve_allocation_subdirectory_components, $staged_basename);
    if ($self->test_name) {
        push @path_components, "test_".$self->test_name;
    }

    push @path_components, $aligner_name_tag;

    my $aligner_version_tag = $self->aligner_version;
    $aligner_version_tag =~ s/[^\w]/_/g;
    push @path_components, $aligner_version_tag;

    if ($self->aligner_params) {
        my $aligner_params_tag = $self->aligner_params;
        $aligner_params_tag =~ s/[^\w]/_/g;
        push @path_components, $aligner_params_tag;
    }

    my $directory = join('/', @path_components);

    $self->debug_message(sprintf("Resolved allocation subdirectory to %s", $directory));
    return $directory;
}

sub _resolve_allocation_subdirectory_components {
    my $self = shift;

    Carp::confess('Class must implement _resolve_allocation_subdirectory_components');
}

1;
