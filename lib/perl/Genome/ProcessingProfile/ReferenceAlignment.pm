package Genome::ProcessingProfile::ReferenceAlignment;

use strict;
use warnings;

use Genome;

class Genome::ProcessingProfile::ReferenceAlignment {
    is => 'Genome::ProcessingProfile::Staged',
    is_abstract => 1,
    subclassify_by => 'subclass_name',
    has => [
        subclass_name => { is_mutable => 0,
                           calculate_from => ['sequencing_platform'],
                           calculate => sub {
                                            my($sequencing_platform) = @_;
                                            Carp::confess "No sequencing platform given to resolve subclass name" unless $sequencing_platform;
                                            return 'Genome::ProcessingProfile::ReferenceAlignment::'.Genome::Utility::Text::string_to_camel_case($sequencing_platform);
                                          }
                         },
    ],

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
        multi_read_fragment_strategy => {
            doc => '',
            is_optional => 1,
        },
        picard_version => {
            doc => 'picard version for MarkDuplicates, MergeSamfiles, CreateSequenceDictionary...',
            is_optional => 1,
        },
        samtools_version => {
            doc => 'samtools version for SamToBam, samtools merge, etc...',
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
        read_calibrator_name => {
            doc => '',
            is_optional => 1,
        },
        read_calibrator_params => {
            doc => '',
            is_optional => 1,
        },
        coverage_stats_params => {
            doc => 'parameters necessary for generating reference coverage in the form of two comma delimited lists split by a colon like 1,5,10,15,20:0,200,500',
            is_optional => 1,
        },
        prior_ref_seq => {
            doc => '',
            is_optional => 1,
        },
        capture_set_name => {
            doc => 'The name of the capture set to evaluate coverage and limit variant calls to within the defined target regions',
            is_optional => 1,
            is_deprecated => 1,
        },
        align_dist_threshold => {
            doc => '',
            is_optional => 1,
        },
    ],
};

sub _resolve_type_name_for_class {
    return 'reference alignment';
}

sub _initialize_build {
    my($self,$build) = @_;

    # Check that the annotator version param is sane before doing the build
    my $annotator_version;
    my $worked = eval {
        my $model = $build->model;
        my $pp = $model->processing_profile;
        $annotator_version = $pp->transcript_variant_annotator_version;
        # When all processing profiles have a param for this, remove this unless block so
        # they'll fail if it's missing
        unless (defined $annotator_version) {
            $annotator_version = Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version;
        }

        my %available_versions = map { $_ => 1 } Genome::Model::Tools::Annotate::TranscriptVariants->available_versions;
        unless ($available_versions{$annotator_version}) {
            die "Requested annotator version ($annotator_version) is not in the list of available versions: "
                . join(', ',keys(%available_versions));
        }
        1;
    };
    unless ($worked) {
        $self->error_message("Could not determine which version of the Transcript Variants annotator to use: $@");
        return;
    }
    return 1;
}

# get alignments (generic name)
sub results_for_instrument_data_input {
    my $self = shift;
    my $input = shift;
    my %segment_info = @_;
    return $self->_fetch_alignment_sets($input,\%segment_info,'get_with_lock');
}

sub results_for_instrument_data_input_with_lock {
    my $self = shift;
    my $input = shift;
    my %segment_info = @_;
    return $self->_fetch_alignment_sets($input,\%segment_info,'get_with_lock');
}

# create alignments (called by Genome::Model::Event::Build::ReferenceAlignment::AlignReads for now...)
sub generate_results_for_instrument_data_input {
    my $self = shift;
    my $input = shift;
    my %segment_info = @_;
    return $self->_fetch_alignment_sets($input,\%segment_info, 'get_or_create');
}

sub _fetch_alignment_sets {
    my $self = shift;
    my $input = shift;
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

sub processing_profile_params_for_alignment {
    my $self = shift;

    my %params = (
                read_aligner_name => $self->read_aligner_name,
                read_aligner_version => $self->read_aligner_version,
                read_aligner_params => $self->read_aligner_params,
                force_fragment => $self->force_fragment,
                read_trimmer_name => $self->read_trimmer_name,
                read_trimmer_version => $self->read_trimmer_version,
                read_trimmer_params => $self->read_trimmer_params,
                picard_version => $self->picard_version,
                samtools_version => $self->samtools_version,
            );

    return \%params;
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
                    filter_name => $input->filter_desc || undef,
                    test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
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
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
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

# TODO: remove
sub prior {
    my $self = shift;
    warn("For now prior has been replaced with the actual column name prior_ref_seq");
    if (@_) {
        die("Method prior() is read-only since it's deprecated");
    }
    return $self->prior_ref_seq();
}

# TODO: remove
sub filter_ruleset_name {
    #TODO: move into the db so it's not constant
    'basic'
}

# TODO: remove
sub filter_ruleset_params {
    ''
}


#< SUBCLASSING >#
#
# This is called by the infrastructure to appropriately classify abstract processing profiles
# according to their type name because of the "sub_classification_method_name" setting
# in the class definiton...
sub _X_resolve_subclass_name {
    my $class = shift;

    my $sequencing_platform;
    if ( ref($_[0]) and $_[0]->can('params') ) {
        my @params = $_[0]->params;
        my @seq_plat_param = grep { $_->name eq 'sequencing_platform' } @params;
        if (scalar(@seq_plat_param) == 1) {
            $sequencing_platform = $seq_plat_param[0]->value;
        }

    }  else {
        my %params = @_;
        $sequencing_platform = $params{sequencing_platform};
    }

    unless ( $sequencing_platform ) {
        my $rule = $class->define_boolexpr(@_);
        $sequencing_platform = $rule->value_for('sequencing_platform');
    }

    return ( defined $sequencing_platform )
    ? $class->_resolve_subclass_name_for_sequencing_platform($sequencing_platform)
    : undef;
}

sub _resolve_subclass_name_for_sequencing_platform {
    my ($class,$sequencing_platform) = @_;
    my @type_parts = split(' ',$sequencing_platform);

    my @sub_parts = map { ucfirst } @type_parts;
    my $subclass = join('',@sub_parts);

    my $class_name = join('::', 'Genome::ProcessingProfile::ReferenceAlignment' , $subclass);
    return $class_name;
}

sub _resolve_sequencing_platform_for_class {
    my $class = shift;

    my ($subclass) = $class =~ /^Genome::ProcessingProfile::ReferenceAlignment::([\w\d]+)$/;
    return unless $subclass;

    return lc join(" ", ($subclass =~ /[a-z\d]+|[A-Z\d](?:[A-Z\d]+|[a-z]*)(?=$|[A-Z\d])/gx));

    my @words = $subclass =~ /[a-z\d]+|[A-Z\d](?:[A-Z\d]+|[a-z]*)(?=$|[A-Z\d])/gx;
    return lc(join(" ", @words));
}

#### IMPLEMENTATION #####

sub stages {
    my $self = shift;
    my $build = shift;
    ## second parameter of each pair is the required flag
    ## if it is 1 and no job events are made at start time
    ## a warning will be printed to the user
    my @stages = (
        reference_preparation   => 1,
        alignment               => 1,
        merge_and_deduplication => 1,
    );

    # Suppress warning for models that do not have reference_coverage_objects,
    # e.g. a region_of_interest_set_name.
    # It will take more work to ensure $build is passed in due to hack used to
    # implement lane QC.
    if (!$build || $self->reference_coverage_objects($build->model)) {
        push @stages, 'reference_coverage' => 1;
    }

    push @stages, (
        variant_detection       => 1,
        transcript_annotation   => 0,
        generate_reports        => 0,
    );

    my @filtered_stages;
    for (my $i=0; $i < $#stages; $i += 2) {
        my $method = $stages[$i] . '_job_classes';

        push @filtered_stages, $stages[$i] if ($stages[$i+1] || $self->$method());
    }

    return @filtered_stages;
}

sub reference_preparation_job_classes {
    my $self = shift;

    my @sub_command_classes;

    my $aligner_class = 'Genome::InstrumentData::AlignmentResult::' . Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($self->read_aligner_name);

    if ($aligner_class->can('prepare_reference_sequence_index')) {
        push @sub_command_classes, 'Genome::Model::Event::Build::ReferenceAlignment::PrepareReferenceSequenceIndex'
    }
    return @sub_command_classes;
}

sub alignment_job_classes {
    my @sub_command_classes= qw/
        Genome::Model::Event::Build::ReferenceAlignment::AlignReads
        Genome::Model::Event::Build::ReferenceAlignment::BamQc
    /;
    return @sub_command_classes;
}

sub reference_coverage_job_classes {
    my $self = shift;
    my $model = shift;
    if ($self->dna_type eq 'cdna' || $self->dna_type eq 'rna') {
        #TODO this needs to be changed to reference build
        my $reference_sequence_build = $model->reference_sequence_build;
        if ($reference_sequence_build->name =~ /^XStrans_adapt_smallRNA_ribo/i) {
            my @steps = (
                'Genome::Model::Event::Build::ReferenceAlignment::RefCov',
            );
            return @steps;
        }
    }
    my @steps = (
        'Genome::Model::Event::Build::ReferenceAlignment::CoverageStats',
    );
    return @steps;
}

sub variant_detection_job_classes {
    my $self = shift;
    my @steps = (
        'Genome::Model::Event::Build::ReferenceAlignment::DetectVariants'
    );
    if(defined $self->snv_detection_strategy || defined $self->indel_detection_strategy ||
            defined $self->sv_detection_strategy || defined $self->cnv_detection_strategy) {
        return @steps;
    }
    else {
        return;
    }
}

sub merge_and_deduplication_job_classes {
    my $self = shift;

    # this is hackish, but maq is a special case right now and hopefully it won't be here long
    if ($self->read_aligner_name eq 'maq') {
        return ('Genome::Model::Event::Build::ReferenceAlignment::DeduplicateLibraries');
    }

    my @steps = (
        'Genome::Model::Event::Build::ReferenceAlignment::MergeAlignments',
    );
    if(defined $self->merger_name) {
        return @steps;
    }
    else {
        return;
    }
}

sub transcript_annotation_job_classes{
    my $self = shift;
    my @steps = (
        'Genome::Model::Event::Build::ReferenceAlignment::AnnotateAdaptor',
        'Genome::Model::Event::Build::ReferenceAlignment::AnnotateTranscriptVariants',
        #'Genome::Model::Event::Build::ReferenceAlignment::AnnotateTranscriptVariantsParallel',
    );
    return @steps;
}

sub generate_reports_job_classes {
    my $self = shift;
    my @steps = (
        'Genome::Model::Event::Build::ReferenceAlignment::RunReports'
    );
    if((defined $self->snv_detection_strategy || defined $self->indel_detection_strategy) && defined $self->duplication_handler_name) {
        return @steps;
    }
    else {
        return;
    }
}

sub reference_preparation_objects {
    "all_sequences";
}

sub alignment_objects {
    my ($self, $model) = @_;

    my @instrument_data = $model->instrument_data;
    my @instrument_data_output = grep {! $_->can('get_segments')} @instrument_data;
    my @segmentable_data = grep {$_->can('get_segments')} @instrument_data;

    for my $instr (@segmentable_data) {
        my @segments = $instr->get_segments();

        # take imported instrument data and split them by read group
        if (@segments > 0 && $self->read_aligner_name ne 'imported' && $instr->isa('Genome::InstrumentData::Imported')) {
            for my $seg (@segments) {
                push @instrument_data_output, {object=>$instr, segment=>$seg};
            }
        } else {
            push @instrument_data_output, $instr;
        }
    }

    return @instrument_data_output;
}

sub reference_coverage_objects {
    my $self = shift;
    my $model = shift;

    my $reference_sequence_build = $model->reference_sequence_build;
    if ($reference_sequence_build->name =~ /^XStrans_adapt_smallRNA_ribo/i) {
        return 'all_sequences';
    }
    my @inputs = Genome::Model::Input->get(model_id => $model->id, name => 'region_of_interest_set_name');
    unless (@inputs) { return; }
    return 'all_sequences';
}


sub variant_detection_objects {
    my $self = shift;
    my $model = shift;
    return 'all_sequences';
}

sub merge_and_deduplication_objects {
    my $self = shift;
    my $model = shift;
    return 'all_sequences';
}

sub generate_reports_objects {
    my $self = shift;
    my $model = shift;
    return 'all_sequences';
}

sub transcript_annotation_objects {
    my $self = shift;
    my $model = shift;
    return unless $model->annotation_reference_build;
    return 'all_sequences';
}

sub default_profile_id {
    return 2635769;
    #This should be updated everytime when new autocron ref-aling default pp is in place
    #The current default is Nov 2011 Default Reference Alignment
}

sub default_profile {
    return __PACKAGE__->get(shift->default_profile_id);
}


1;
