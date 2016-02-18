package Genome::ProcessingProfile::ReferenceAlignment;

use strict;
use warnings;

use Genome;

class Genome::ProcessingProfile::ReferenceAlignment {
    is => 'Genome::ProcessingProfile',
    has_param => [
        sequencing_platform => {
            doc => 'The sequencing platform from whence the model data was generated',
            valid_values => ['454', 'solexa', 'sanger'],
        },
        dna_type => {
            doc => 'the type of dna used in the reads for this model',
            valid_values => ['genomic dna', 'cdna']
        },
        transcript_variant_annotator_version => {
            doc => 'Version of the "annotate transcript-variants" tool to run during the annotation step',
            is_optional => 1,
            default_value => Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version,
            valid_values => [ 0,1,2,3,4],#Genome::Model::Tools::Annotate::TranscriptVariants->available_versions ],
        },
        transcript_variant_annotator_filter => {
            doc => 'annotation-filter option to be used by the "annotate transcript-variants" tool run during the annotation step',
            is_optional => 1,
            default_value => 'top',
            valid_values => ['top', 'none', 'gene'],
        },
        transcript_variant_annotator_accept_reference_IUB_codes => {
            doc => 'annotation accept-reference-IUB-codes option to be used by the "annotate transcript-variants" to run during the annotation step',
            is_optional => 1,
            default_value => '0',
            valid_values => [0, 1],
        },
        alignment_strategy => {
            is => 'Text',
            is_many => 0,
            is_optional => 1,
            doc => 'Strategy to be used to align the instrument data (this will eventually replace many of the alignment parameters)',
        },
        snv_detection_strategy => {
            is => "Text",
            is_many => 0,
            is_optional =>1,
            doc => "Strategy to be used to detect snvs.",
        },
        indel_detection_strategy => {
            is => "Text",
            is_many => 0,
            is_optional =>1,
            doc => "Strategy to be used to detect indels.",
        },
        sv_detection_strategy => {
            is => "Text",
            is_many => 0,
            is_optional =>1,
            doc => "Strategy to be used to detect svs.",
        },
        cnv_detection_strategy => {
            is => "Text",
            is_many => 0,
            is_optional =>1,
            doc => "Strategy to be used to detect cnvs.",
        },
        picard_version => {
            doc => 'picard version for MarkDuplicates, MergeSamfiles, CreateSequenceDictionary...',
            is_optional => 1,
        },
        samtools_version => {
            doc => 'samtools version for SamToBam, samtools merge, etc...',
            is_optional => 1,
        },
        bedtools_version => {
            doc => 'bedtools version for bedtools bamtofastq',
            is_optional => 1,
        },
        merger_name => {
            doc => 'name of bam merger, picard, samtools (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        merger_version => {
            doc => 'version of bam merger (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        merger_params => {
            doc => 'parameters of bam merger (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        duplication_handler_name => {
            doc => 'name of pcr duplication handler (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        duplication_handler_version => {
            doc => 'version of pcr duplication handler (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        duplication_handler_params => {
            doc => 'parameters of pcr duplication handler (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        read_aligner_name => {
            doc => 'alignment algorithm/software used for this model (this will be replaced by alignment_strategy)',
        },
        read_aligner_version => {
            doc => 'the aligner version used for this model (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        read_aligner_params => {
            doc => 'command line args for the aligner (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        force_fragment => {
            is => 'Integer',
            #This doesn't seem to work yet because of the create code, can't the valid values logic be removed from create???
            default_value => '0',
            #valid_values => ['0', '1'],
            doc => 'force all alignments as fragment reads',
            is_optional => 1,
        },
        read_trimmer_name => {
            doc => 'trimmer algorithm/software used for this model',
            is_optional => 1,
        },
        read_trimmer_version => {
            doc => 'the trimmer version used for this model',
            is_optional => 1,
        },
        read_trimmer_params => {
            doc => 'command line args for the trimmer',
            is_optional => 1,
        },
        coverage_stats_params => {
            doc => 'parameters necessary for generating reference coverage in the form of two comma delimited lists split by a colon like 1,5,10,15,20:0,200,500',
            is_optional => 1,
        },
        append_event_steps => {
            doc => 'Event classes to append to event_stage_job_classes, e.g. "alignment => Genome::Model::Event::Build::ReferenceAlignment::QC::CopyNumber".',
            is_optional => 1,
            is_deprecated => 1,
        },

    ],
};

#TODO Once old snapshots stop subclassifying by sequencing platform, remove this!
sub get {
    my $class = shift;

    if(ref $class) {
        return $class->SUPER::get(@_);
    }

    my $bx = UR::BoolExpr->resolve_normalized('Genome::ProcessingProfile::ReferenceAlignment', @_);
    #until the database is updated, need to check all possible class names
    $bx = $bx->add_filter(subclass_name => ['Genome::ProcessingProfile::ReferenceAlignment', 'Genome::ProcessingProfile::ReferenceAlignment::Solexa', 'Genome::ProcessingProfile::ReferenceAlignment::454'],);

    return Genome::ProcessingProfile->get($bx);
}

sub _resolve_type_name_for_class {
    return 'reference alignment';
}

# get alignments (generic name)
sub results_for_instrument_data_input {
    my $self = shift;
    my $input = shift;
    my $result_users = shift;
    my %segment_info = @_;
    return $self->_fetch_alignment_sets($input,$result_users,\%segment_info,'get_with_lock');
}

sub results_for_instrument_data_input_with_lock {
    my $self = shift;
    my $input = shift;
    my $result_users = shift;
    my %segment_info = @_;
    return $self->_fetch_alignment_sets($input,$result_users,\%segment_info,'get_with_lock');
}

# create alignments (called by Genome::Model::Event::Build::ReferenceAlignment::AlignReads for now...)
sub generate_results_for_instrument_data_input {
    my $self = shift;
    my $input = shift;
    my $result_users = shift;
    my %segment_info = @_;
    return $self->_fetch_alignment_sets($input,$result_users,\%segment_info, 'get_or_create');
}

sub _fetch_alignment_sets {
    my $self = shift;
    my $input = shift;
    my $result_users = shift;
    my $segment_info = shift;
    my $mode = shift;

    my $model = $input->model;

    my @param_sets = $self->params_for_alignment($input);
    unless (@param_sets) {
        $self->error_message('Could not get alignment parameters for this instrument data input');
        return;
    }
    my @alignments;
    for (@param_sets)  {
        my %params = %$_;
        $params{users} = $result_users;

        # override segments if requested
        if (exists $segment_info->{instrument_data_segment_id}) {
            delete $params{instrument_data_segment_id};
            delete $params{instrument_data_segment_type};
        }
        my $alignment = Genome::InstrumentData::AlignmentResult->$mode(%params, %$segment_info);
        unless ($alignment) {
             #$self->error_message("Failed to $mode an alignment object");
             return;
         }
        push @alignments, $alignment;
    }
    return @alignments;
}

sub params_for_alignment {
    my $self = shift;
    my $input = shift;

    my $model = $input->model;
    my $reference_build = $model->reference_sequence_build;
    my $reference_build_id = $reference_build->id;

    unless ($self->type_name eq 'reference alignment') {
        $self->error_message('Can not create an alignment object for model type '. $self->type_name);
        return;
    }

    my %params = (
                    instrument_data_id => $input->value_id || undef,
                    aligner_name => $self->read_aligner_name || undef,
                    reference_build_id => $reference_build_id || undef,
                    aligner_version => $self->read_aligner_version || undef,
                    aligner_params => $self->read_aligner_params || undef,
                    force_fragment => $self->force_fragment || undef,
                    trimmer_name => $self->read_trimmer_name || undef,
                    trimmer_version => $self->read_trimmer_version || undef,
                    trimmer_params => $self->read_trimmer_params || undef,
                    picard_version => $self->picard_version || undef,
                    samtools_version => $self->samtools_version || undef,
                    bedtools_version => $self->bedtools_version || undef,
                    filter_name => $input->filter_desc || undef,
                    test_name => Genome::Config::get('software_result_test_name') || undef,
                    instrument_data_segment_type => undef,
                    instrument_data_segment_id => undef,
                );

    my @param_set = (\%params);
    return @param_set;
}

sub params_for_merged_alignment {
    my $self = shift;
    my $build = shift; #TODO possibly calculate segment info in _fetch_merged_alignment_result (which calls this)
    my @inputs = @_;

    my $filters = [];
    for my $i (0..$#inputs) {
        my $input = $inputs[$i];
        if($input->filter_desc) {
            push @$filters, join(':', $input->value->id, $input->filter_desc);
        }
    }

    my $segment_parameters = [];
    if($build) {
        my @align_reads_events = grep {$_->isa('Genome::Model::Event::Build::ReferenceAlignment::AlignReads')} $build->events;
        for my $i (0..$#inputs) {
            my $input = $inputs[$i];
            my @alignment_events = grep {$_->instrument_data_id eq $input->value->id} @align_reads_events;

            #if multiple events, this is a chunked alignment
            if (@alignment_events > 1 or (@alignment_events == 1 and defined $alignment_events[0]->instrument_data_segment_id)) {
                for my $alignment_event (@alignment_events) {
                    push @$segment_parameters, join(':', $alignment_event->instrument_data_id, $alignment_event->instrument_data_segment_id, $alignment_event->instrument_data_segment_type);
                }
            }
        }
    }

    my $instrument_data = [];

    for my $i (0..$#inputs) {
        push @$instrument_data, $inputs[$i]->value->id;
    }

    my %params = (
        instrument_data_id => $instrument_data,
        reference_build_id => $build->reference_sequence_build->id,
        merger_name => $self->merger_name,
        merger_params => $self->merger_params,
        merger_version => $self->merger_version,
        duplication_handler_name => $self->duplication_handler_name,
        duplication_handler_params => $self->duplication_handler_params,
        duplication_handler_version => $self->duplication_handler_version,
        aligner_name => $self->read_aligner_name || undef,
        aligner_version => $self->read_aligner_version || undef,
        aligner_params => $self->read_aligner_params || undef,
        force_fragment => $self->force_fragment || undef,
        trimmer_name => $self->read_trimmer_name || undef,
        trimmer_version => $self->read_trimmer_version || undef,
        trimmer_params => $self->read_trimmer_params || undef,
        picard_version => $self->picard_version || undef,
        samtools_version => $self->samtools_version || undef,
        bedtools_version => $self->bedtools_version || undef,
        test_name => Genome::Config::get('software_result_test_name') || undef,
    );
    if(scalar @$filters) {
        $params{filter_name} = $filters;
    }
    if(scalar @$segment_parameters) {
        $params{instrument_data_segment} = $segment_parameters;
    }

    my @param_set = (\%params);
    return @param_set;
}

1;
