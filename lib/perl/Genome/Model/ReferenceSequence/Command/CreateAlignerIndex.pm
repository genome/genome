package Genome::Model::ReferenceSequence::Command::CreateAlignerIndex;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceSequence::Command::CreateAlignerIndex {
    is => ['Command::V2'],
    has_input => [
       aligner_name => {
           is => "Text",
       },
       aligner_version => {
           is => "Text",
       },
       aligner_params => {
           is => "Text",
       },
       reference_sequence_build_id => {
           is => 'Text',
       },
       annotation_build_id => {
           is => 'Text',
           is_optional => 1,
       },
    ],
    has => [
       reference_sequence_build => {
           id_by => 'reference_sequence_build_id',
           is => 'Genome::Model::Build::ReferenceSequence',
       },
       annotation_build => {
            id_by => 'annotation_build_id',
            is => 'Genome::Model::Build::ImportedAnnotation',
            is_optional => 1,
       },
    ],
    has_param => [
        lsf_resource => {
            is => 'Text',
            value => "-R 'select[type==LINUX64 && mem>16000 && tmp>150000] span[hosts=1] rusage[tmp=150000, mem=16000]' -M 16000000",
        }
    ],
};

sub alignment_result_class {
    my $self = shift;
    return 'Genome::InstrumentData::AlignmentResult::' . Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($self->aligner_name);
}

sub bsub_rusage {
    my $self = shift;
    my $delegate = $self->alignment_result_class;
    my $rusage = $delegate->required_rusage_for_building_index(reference_build=>$self->reference_sequence_build);
    return $rusage;
}

sub execute {
    my $self = shift;

    return $self->_process('get_or_create');
}


sub shortcut {
    my $self = shift;

    return $self->_process('get_with_lock');
}

sub _process {
    my $self = shift;
    my $mode = shift;

    my @errors;

    my %params_for_reference = (
        aligner_name=>$self->aligner_name,
        aligner_params=>$self->aligner_params,
        aligner_version=>$self->aligner_version,
        reference_build=>$self->reference_sequence_build,
    );

    $self->debug_message(sprintf("Finding or generating reference build index for aligner %s version %s params %s refbuild %s ",
                                                $self->aligner_name, $self->aligner_version,
                                                $self->aligner_params, $self->reference_sequence_build->id));

    my $reference_sequence_index;
    my $annotation_index;
    if ($mode eq 'get_or_create') {
        $reference_sequence_index = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_or_create(%params_for_reference);
        unless ($reference_sequence_index) {
            $self->error_message("Error getting or creating reference sequence index!");
            push @errors, $self->error_message;
        }

        if($self->annotation_build_id) {
            $annotation_index = Genome::Model::Build::ReferenceSequence::AnnotationIndex->get_or_create(%params_for_reference, annotation_build => $self->annotation_build);
            unless($annotation_index) {
                $self->error_message('Error getting or creating annotation index!');
                push @errors, $self->error_message;
            }
        }
    } elsif ($mode eq 'get_with_lock') {
        $reference_sequence_index = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_with_lock(%params_for_reference);
        unless ($reference_sequence_index) {
            return undef;
        }

        if($self->annotation_build_id) {
            $annotation_index = Genome::Model::Build::ReferenceSequence::AnnotationIndex->get_with_lock(%params_for_reference, annotation_build => $self->annotation_build);
            unless($annotation_index) {
                return undef;
            }
        }
    } else {
        $self->error_message("generate alignment reference index mode unknown: $mode");
        die $self->error_message;
    }

    if (@errors) {
        $self->error_message(join("\n",@errors));
        if($mode eq 'get') {
            return 0;
        } else {
            die $self->error_message;
        }
    }

    $self->debug_message("Complete!");
    return 1;
}

1;

