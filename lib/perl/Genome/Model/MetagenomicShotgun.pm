package Genome::Model::MetagenomicShotgun;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun {
    is => 'Genome::ModelDeprecated',
    has_param => [
        filter_contaminant_fragments => {
            is => 'Boolean',
            doc => 'when set, reads with mate mapping to contamination reference are considered contaminated, and not passed on to subsequent alignments',
            default => 0,
        },
        filter_duplicates => {
            is => 'Boolean',
            doc =>  'when set, duplicate reads are filtered out when extracting unaligned reads from contamination screen alignment',
        },
        contamination_screen_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for contamination screen',
            is_optional=> 1,
        },
        metagenomic_nucleotide_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for metagenomic alignment',
        },
        metagenomic_protein_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for realignment of unaligned reads from first metagenomic alignment',
            is_optional => 1,
        },
        viral_nucleotide_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for first viral verification alignment',
            is_optional => 1,
        },
        viral_protein_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for first viral verification alignment',
            is_optional => 1,
        },
    ],
    has => [
        contamination_screen_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            is_optional => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'contamination_screen_reference', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        metagenomic_protein_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'metagenomic_protein_reference ', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        viral_nucleotide_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            is_optional => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'viral_nucleotide_reference ', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        viral_protein_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            is_optional => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'viral_protein_reference ', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        metagenomic_nucleotide_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'metagenomic_alignment_reference', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        contamination_screen_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1,
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'contamination_screen_model'],
        },
        metagenomic_nucleotide_model => {
            is => 'Genome::Model::ReferenceAlignment',
            via => 'from_model_links', 
            to => 'from_model',
            where => [role => 'metagenomic_nucleotide_model'],
        },
        metagenomic_protein_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1,
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'metagenomic_protein_model'],
        },
        viral_nucleotide_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1, 
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'viral_nucleotide_model'],
        },
        viral_protein_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1, 
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'viral_protein_model'],
        },
    ],
};

sub sub_model_labels {
    return (qw/ contamination_screen metagenomic_nucleotide metagenomic_protein viral_nucleotide viral_protein /);
}

sub build_subclass_name {
    return 'metagenomic-composition-shotgun';
}

sub delete {
    my $self = shift;
    for my $sub_model ($self->from_models) {
        $sub_model->delete;
    }
    return $self->SUPER::delete(@_);
}

sub create {
    my ($class, %params) = @_;

    my $self = $class->SUPER::create(%params);
    return unless $self;

    my $processing_profile = $self->processing_profile;
    for ( $self->sub_model_labels ){
        my $pp_method = $_."_pp_id";
        if($self->processing_profile->$pp_method) {
            my $model = $self->_create_model_for_type($_);
            unless ($model) {
                $self->error_message("Error creating $_ model!");
                $self->delete;
                return;
            }
        }
    }

    return $self;
}

sub sequencing_platform{
    return 'solexa';
}

# TODO make this get or create
sub _create_model_for_type {
    my $self = shift;
    my $type = shift;

    #CREATE UNDERLYING REFERENCE ALIGNMENT MODELS
    my $pp_accessor = $type."_pp_id";
    my $reference_accessor = $type."_reference";
    my %model_params = (
        processing_profile_id => $self->processing_profile->$pp_accessor,
        subject_name => $self->subject_name,
        name => $self->name.".$type model",
        reference_sequence_build=>$self->$reference_accessor,
    );
    my $model = Genome::Model::ReferenceAlignment->create( %model_params );

    unless ($model){
        die $self->error_message("Couldn't create contamination screen model with params ".join(", ", map {$_ ."=>". $model_params{$_}} keys %model_params) );
    }

    $self->add_from_model(from_model=> $model, role=>$type.'_model');
    $self->status_message("Created $type model ".$model->__display_name__);

    return $model;
}

sub _resolve_resource_requirements_for_build {
    my ($self, $build) = @_;
    my @instrument_data = $build->instrument_data;
    my $gtmp = 30 + 5 * (1 + scalar(@instrument_data));
    return "-R 'rusage[gtmp=$gtmp:mem=16000]' -M 16000000";
}

sub sub_model_for_label {
    my ($self, $sub_model_label) = @_;

    if ( not $sub_model_label ) {
        Carp::confess( $self->error_message('No sub model label given!') ) 
    }

    if ( not grep { $sub_model_label eq $_ } $self->sub_model_labels ) {
        Carp::confess( $self->error_message('Invalid sub model label! '.$sub_model_label) );
    }

    my $sub_model_method = $sub_model_label.'_model';
    my $sub_model = $self->$sub_model_method;
    if ( not $sub_model ) {
        Carp::confess( $self->error_message("Failed to get sub model ($sub_model_label) for meta shot model ".$self->__display_name__) );
    }

    return $sub_model;
}

sub last_complete_build_for_sub_model {
    my ($self, $sub_model_label) = @_;

    my $sub_model = $self->sub_model_for_label($sub_model_label); # confess
    my $sub_build = $sub_model->last_complete_build;
    if ( not $sub_build ) {
        $self->error_message("Failed to get sub build for sub model ($sub_model_label) for meta shot model ".$self->__display_name__);
        return;
    }

    return $sub_build;
}

#< Work Flow >#
sub map_workflow_inputs {
    my ($self, $build) = @_;
    my @instrument_data = $build->instrument_data;
    return (
        input_build => $build,
        instrument_data => \@instrument_data,
        # contamination screen
        align_to_contamination_screen => 'contamination_screen',
        extract_from_contamination_screen => 'unaligned',
        # meta nt
        align_to_meta_nt => 'metagenomic_nucleotide',
        extract_from_meta_nt => [qw/ aligned unaligned /],
        # meta nr
        align_to_meta_nr => 'metagenomic_protein',
        extract_from_meta_nr => 'aligned',
        # viral
        align_to_viral => [qw/ viral_nucleotide viral_protein /],
    );
}

sub _resolve_workflow_for_build {
    my ($self, $build, $lsf_queue, $lsf_project) = @_;

    # Align original instrument data to contamination screen reference
    # Extract unaligned reads from contamination screen alignment
    # Align unaligned reads from contamination screen alignment to meta nt reference
    # Extract aligned & unaligned reads from meta nt alignment
    # Align unaligned reads from meta nt alignment to meta nr reference
    # Extract aligned reads from meta nr
    # Align aligned reads from meta nt and meta nr alignments to viral nt and nr references
    # Link Alignments
    # TODO
    # Reports

    # Create work flow and set queue and project
    my $workflow = Workflow::Model->create(
        name => $build->workflow_name,
        input_properties => [qw/ 
            build instrument_data 
            align_to_contamination_screen extract_from_contamination_screen 
            align_to_meta_nt extract_from_meta_nt
            align_to_meta_nr extract_from_meta_nr
            align_to_viral
        /],
        output_properties => [qw/ build /],
        log_dir => $build->log_directory,
    );
    $lsf_queue //= $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT};
    $lsf_project //= 'build' . $build->id;
    my $left_operation = $workflow->get_input_connector;

    # Generic add operation anon sub
    my $add_operation = sub {
        my (%params) = @_;

        my $operation_name = delete $params{operation_name};
        Carp::confess( $self->error_message(
                'Required param "operation_name" missing! Cannot add operation to work flow!'
            ) ) if not $operation_name;
        my $sub_model_type = delete $params{sub_model_type};
        Carp::confess( $self->error_message(
                'Missing param "sub_model_type" to add align to operation to work flow!'
            ) ) if not $sub_model_type;
        my $links_from_left_op = delete $params{links_from_left_op};
        $links_from_left_op ||= {};
        $links_from_left_op->{build} = 'input_build';
        my $links_from_input_connector = delete $params{links_from_input_connector};
        $links_from_input_connector ||= {};
        Carp::confess( $self->error_message(
                'Unknown params to add operation to work flow! '.Data::Dumper::Dumper(\%params)
            ) ) if %params;

        my $command_sub_name = join('', map { ucfirst } split(' ', $operation_name));
        my $command_class_name = 'Genome::Model::MetagenomicShotgun::Build::'.$command_sub_name;
        my $operation_type = Workflow::OperationType::Command->create(command_class_name => $command_class_name);
        Carp::confess( $self->error_message("Failed to create work flow operation for $command_class_name") ) if not $operation_type;
        
        $operation_type->lsf_queue($lsf_queue);
        $operation_type->lsf_project($lsf_project);

        my $operation = $workflow->add_operation(
            name => $operation_name.' '.join(' ', split('_', $sub_model_type)),
            operation_type => $operation_type,
        );

        for my $left_property ( keys %$links_from_left_op ) {
            $workflow->add_link(
                left_operation => $left_operation,
                left_property => $left_property,
                right_operation => $operation,
                right_property => $links_from_left_op->{$left_property},
            );
        }

        my $input_connector = $workflow->get_input_connector;
        for my $left_property ( keys %$links_from_input_connector ) {
            $workflow->add_link(
                left_operation => $input_connector,
                left_property => $left_property,
                right_operation => $operation,
                right_property => $links_from_input_connector->{$left_property},
            );
        }

        return $operation;
    };

    # Align to operation builder
    my $add_align_to_operation = sub {
        my ($sub_model_type) = @_;

        Carp::confess( $self->error_message('Missing param "sub_model_type" to add align to operation to work flow!') ) if not $sub_model_type;

        my $operation = $add_operation->(
            operation_name => 'align to',
            sub_model_type => $sub_model_type,
            links_from_left_op => { instrument_data => 'instrument_data', },
            links_from_input_connector => { 'align_to_'.$sub_model_type => 'sub_model_label', },
        ); # confesses
        $operation->parallel_by('sub_model_label') if $sub_model_type eq 'viral';

        return $operation;
    };

    # Extract from alignment operation builder
    my $add_extract_from_operation = sub {
        my ($sub_model_type) = @_;

        Carp::confess( $self->error_message('Missing param "sub_model_type" to add align to operation to work flow!') ) if not $sub_model_type;

        my $operation = $add_operation->(
            operation_name => 'extract from',
            sub_model_type => $sub_model_type,
            links_from_left_op => {
                'sub_model_label' => 'sub_model_label',
            },
            links_from_input_connector => {
                'extract_from_'.$sub_model_type => 'type',
            },
        ); # confesses
        $operation->parallel_by('type') if $sub_model_type =~ /^meta/;

        return $operation;
    };

    # Align original instrument data to contamination screen reference
    my $align_to_contamination_screen_reference_op = $add_align_to_operation->('contamination_screen'); #confess
    $left_operation = $align_to_contamination_screen_reference_op;

    # Extract unaligned reads from contamination screen alignment
    my $extract_unaligned_from_contamination_screen_op = $add_extract_from_operation->('contamination_screen'); #confess
    $left_operation = $extract_unaligned_from_contamination_screen_op;

    # Align unaligned reads from contamination screen alignment to meta nt reference
    my $align_to_meta_nt_reference_op = $add_align_to_operation->('meta_nt'); #confess
    $left_operation = $align_to_meta_nt_reference_op;

    # Extract aligned & unaligned reads from meta nt alignment
    my $extract_aligned_and_unaligned_from_meta_nt_op = $add_extract_from_operation->('meta_nt');
    $left_operation = $extract_aligned_and_unaligned_from_meta_nt_op;

    # Align unaligned reads from meta nt alignment to meta nr reference
    my $align_to_meta_nr_op = $add_align_to_operation->('meta_nr');
    $left_operation = $align_to_meta_nr_op;

    # Extract aligned & unaligned reads from meta nr
    my $extract_aligned_and_unaligned_from_meta_nr_op = $add_extract_from_operation->('meta_nr'); #confess
    $left_operation = $extract_aligned_and_unaligned_from_meta_nr_op;

    # Align aligned reads from meta nt and meta nr alignments to viral nt and nr references
    my $align_to_viral_op = $add_align_to_operation->('viral');
    $left_operation = $align_to_viral_op;
    
    # Link Alignments
    my $link_alignments_op = $add_operation->(operation_name => 'link alignments', sub_model_type => 'all');
    $left_operation = $link_alignments_op;

    # Output link
    $workflow->add_link(
        left_operation => $left_operation,
        left_property => 'build',
        right_operation => $workflow->get_output_connector,
        right_property => 'build',
    );

    return $workflow;
}
#<>#

1;

