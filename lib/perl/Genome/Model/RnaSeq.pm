package Genome::Model::RnaSeq;

use strict;
use warnings;

use Genome;
use version;
use Genome::Utility::Text;
use File::Spec;

class Genome::Model::RnaSeq {
    is => 'Genome::ModelDeprecated',
    has => [
        subject                      => { is => 'Genome::Sample', id_by => 'subject_id' },
        reference_sequence_name      => { via => 'reference_sequence_build', to => 'name' },
    ],
    has_input => [
        reference_sequence_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
        },
        annotation_build => {
            is => "Genome::Model::Build::ImportedAnnotation",
            is_optional => 1,
        },
        target_region_set_name => {
            is => 'Text',
            is_optional => 1,
            doc => 'limits the assignment of instrument data by default to only data with a matching TRSN'
        },
        cancer_annotation_db => {
            is => 'Text',
            is_optional => 1,
            doc => 'db of cancer annotation (see \'genome db list\' for latest version of desired database)',
        },
    ],
    has_param => [
        sequencing_platform => {
            doc => 'The sequencing platform from whence the model data was generated',
            valid_values => ['454', 'solexa'],
            is_optional => 1,
        },
        digital_expression_detection_strategy => {
            is => 'Text',
            is_optional => 1,
            example_values => ['htseq-count 0.5.4p1 [--mode intersect-strict --minaqual 1 --blacklist-alignments-flags 0x0104 --results-version 1]',],
            doc => 'measure expression level using exact read counts (non-normalized)'
        },
        dna_type => {
            # TODO: most of the samples are now flagged as 'rna' ...so this field probably isn't doing anything -ssmith
            doc => 'the type of dna used in the reads for this model',
            valid_values => ['cdna'],
            is_optional => 1,
        },
        read_aligner_name => {
            doc => 'alignment algorithm/software used for this model',
            is_optional => 1,
        },
        read_aligner_version => {
            doc => 'the aligner version used for this model',
            is_optional => 1,
        },
        read_aligner_params => {
            doc => 'command line args for the aligner',
            is_optional => 1,
        },
        expression_name => {
            doc => 'algorithm used to detect expression levels',
            is_optional => 1,
        },
        expression_version => {
            doc => 'the expression detection version used for this model',
            is_optional => 1,
        },
        expression_params => {
            doc => 'the expression detection params used for this model',
            is_optional => 1,
        },
        samtools_version => {
            doc => 'the version of Samtools to use when manipulating SAM/BAM files',
            is_optional => 1,
        },
        picard_version => {
            doc => 'the version of Picard to use when manipulating SAM/BAM files',
            is_optional => 1,
        },
        picard_strand_specificity => {
            doc => 'The transcript strand used by Picard for evaluatin RNA-seq QC metrics.',
            valid_values => Genome::Model::Tools::Picard::CollectRnaSeqMetrics->__meta__->property("strand_specificity")->valid_values,
            is_optional => 1,
        },
        deduplication_handler => {
            doc => 'the software used to deduplicate or mark duplicates in the aligned BAM file',
            is_optional => 1,
        },
        samstat_version => {
            doc => 'the version of Samstat to use for BamQc',
            is_optional => 1,
        },
        fastqc_version => {
            doc => 'the version of FastQC to use for BamQc',
            is_optional => 1,
        },
        bedtools_version => {
            doc => 'the version of BEDTools to use for splice junction reporting.',
            is_optional => 1,
        },
        error_rate_version  => {
            doc => 'the version of error rate C tool to use for BamQc.',
            is_optional => 1,
        },
        calculate_error_rate => {
            doc => 'A flag to calculate error rate during BamQc',
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
        annotation_reference_transcripts_mode => {
            doc => 'The mode to use annotation_reference_transcripts for expression analysis',
            is => 'Text',
            is_optional => 1,
        },
        mask_reference_transcripts => {
            doc => 'The mask level to ignore transcripts located in these annotation features',
            is_optional => 1,
            valid_values => ['rRNA','MT','pseudogene','rRNA_MT','rRNA_MT_pseudogene'],
        },
        transcriptome_coverage_annotation_file_basenames => {
            is_optional => 1,
            doc => 'A comma delimited list of annotation file basenames to generate transcriptome coverage metrics for.',
        },
        transcript_coverage_merge_annotation_features => {
            is_optional => 1,
            doc => 'A yes, no or both answer for merging exons for transcriptome coverage metrics.',
        },
        transcriptome_coverage_mask_reference_transcripts => {
            doc => 'The mask level to ignore transcripts located in these annotation features',
            is_optional => 1,
            valid_values => ['rRNA','MT','pseudogene','rRNA_MT','rRNA_MT_pseudogene'],
        },
        fusion_detector => {
            is_optional => 1,
            is => 'Text',
            valid_values => ['chimerascan', 'chimerascan-vrl'],
            doc => 'The program to use for detecting fusion events',
        },
        fusion_detector_version => {
            is_optional => 1,
            is => 'Text',
            doc => 'The version of the fusion-detector to use.',
        },
        fusion_detector_params => {
            is_optional => 1,
            doc => 'Detector specific fusion-detection parameters.',
        },
        bowtie_version => {
            is_optional => 1,
            is => 'Text',
            doc => 'version of bowtie for tophat to use internally',
        },
    ],
    doc => 'A genome model produced by aligning cDNA reads to a reference sequence.',
};

sub compatible_instrument_data {
    my $self = shift;
    my @compatible_instrument_data = $self->SUPER::compatible_instrument_data(@_);
    return grep{!($_->can('is_paired_end')) or $_->is_paired_end} @compatible_instrument_data;
}

sub _parse_strategy {
    my ($self,$strategy,$build,$label) = @_;

    # TODO: replace this with a real parser, and pull up into the general model class.
    # This will not handle boolean logic or many odd characters
    # It will let us make forward-compatible specs for tool/version/params much like DV2 for now.

    my ($tool, $version, $x, $params);
    unless ( ($tool, $version, $x, $params) = ($strategy =~ /^(.+?)\s+(\S+)\s*(|\[(.*)\])\s*$/) ) {
        die "failed to parse strategy: $strategy!";
    }

    #Genome::Utility::Text::string_to_camel_case($detector,"-")
    my $cmd_class = 'Genome::Model::Tools::' . join('::', map { ucfirst(lc($_)) } split('-',$tool));

    my %params;
    my %params1 = ($params ? split(/\s+/,$params) : ());
    for my $param (keys %params1) {
        my $value = $params1{$param};
        unless ($param =~ s/^--//) {
            die "unexpected param $param in strategy: expected -- prefix";
        }
        $param =~ s/-/_/g;
        $params{$param} = $value;
    }

    if ($cmd_class->can('app_version')) {
        $params{app_version} = $version;
    }
    elsif ($cmd_class->can('use_version')) {
        $params{use_version} = $version;
    }

    return $cmd_class, \%params;
}

sub map_workflow_inputs {
    my ($self, $build) = @_;

    Carp::confess('No build given to map workflow inputs!') if not $build;

    # modes are validated later in Genome::Model::Build::RnaSeq and
    # Genome::InstrumentData::AlignmentResult::Command::CufflinksExpression
    my @modes = split /\s*,\s*/, $self->processing_profile->annotation_reference_transcripts_mode;

    # remove any extra whitespace
    for (@modes) {s/\s+/ /g; s/^\s+//; s/\s+$//}

    my @inputs = (
        ($build ? (build_id => $build->id) : ()),
        annotation_reference_transcripts_mode => \@modes,
    );

    # add any new strategies to the list below, or make this more fully dynamic
    for my $strategy_name (qw/digital_expression/) {
        my $strategy_value_method = $strategy_name . '_detection_strategy';
        if (my $strategy = $self->$strategy_value_method) {
            my $strategy_parser_method = '_parse_' . $strategy_value_method;
            unless ($self->can($strategy_parser_method)) {
                # no _parse_BLAH_detection_strategy method
                # fall back to the default
                $strategy_parser_method = '_parse_strategy';
            }
            my ($class,$params) = $self->$strategy_parser_method($strategy, $build);

            my $result_users = Genome::SoftwareResult::User->user_hash_for_build($build);

            $params->{sponsor} = $result_users->{sponsor};
            $params->{requestor} = $result_users->{requestor};
            $params->{user} = $build;
            $params->{label} = $strategy_name . '_result';
            $params->{output_dir} = File::Spec->join($build->data_directory, 'results', $strategy_name . '_result');
            for my $key (keys %$params) {
                my $input_name = $strategy_name . '_' . $key;
                my $value = $params->{$key};
                push @inputs, $input_name => $value;
            }
        }
    }

    if ($self->fusion_detector) {
        push @inputs, $self->fusion_detection_inputs($build->processing_profile);
        push @inputs, $self->chimerascan_annotation_inputs($build) if $self->cancer_annotation_db;
    }

    my %inputs = @inputs;

    return @inputs;
}

sub _resolve_workflow_for_build {
    # This is called by Genome::Model::Build::start()
    my $self = shift;
    my $build = shift;
    my $lsf_queue = shift; # TODO: the workflow shouldn't need this yet
    my $lsf_project = shift;

    if (!defined $lsf_queue || $lsf_queue eq '' || $lsf_queue eq 'inline') {
        $lsf_queue = Genome::Config::get('lsf_queue_build_worker_alt');
    }
    if (!defined $lsf_project || $lsf_project eq '') {
        $lsf_project = 'build' . $build->id;
    }

    my ($version_number) = $self->read_aligner_version =~ /^([\d+])/;
    my $aligner_name = $self->read_aligner_name;

    my $processing_profile = $self->processing_profile;

    my $run_splice_junction_summary;

    #SpliceJunctionSummary need junctions.bed produced by tophat aligner
    if ($aligner_name eq 'tophat' && $self->bedtools_version && $self->annotation_build) {
        my $annotation_build = $self->annotation_build;
        my ($ensembl_version) = split(/_/,$annotation_build->version);
        unless ($ensembl_version) {
            die('Failed to parse Ensembl version from annotation build: '. $annotation_build->id);
        }
        if ($ensembl_version >= 67) {
            $run_splice_junction_summary = 1;
        }
        else {
            $self->debug_message('Skipping SpliceJunctionSummary for annotation build: '. $self->annotation_build);
        }
    }

    my $workflow = Genome::WorkflowBuilder::DAG->create(
        name => $build->workflow_name,
    );

    # Alignment
    my %params = ($version_number < 2)? (
        name => 'RnaSeq Tophat Alignment',
        command => 'Genome::Model::RnaSeq::Command::AlignReads::Tophat',
    ) : (
        name => 'RnaSeq Alignment',
        command => 'Genome::Model::RnaSeq::Command::AlignReads',
    );

    my $alignment_operation = Genome::WorkflowBuilder::Command->create(%params);
    $workflow->add_operation($alignment_operation);

    $alignment_operation->lsf_queue($lsf_queue);
    $alignment_operation->lsf_project($lsf_project);

    $workflow->connect_input(
        input_property => 'build_id',
        destination => $alignment_operation,
        destination_property => 'build_id'
    );

    # Digital expression
    my $digital_expression_detection_operation = undef;
    if (my $strategy = $processing_profile->digital_expression_detection_strategy) {
        my ($class,$params) = $self->_parse_strategy($strategy, $build);

        $digital_expression_detection_operation = Genome::WorkflowBuilder::Command->create(
            name => 'RnaSeq Digital Expression Detection',
            command => $class,
        );
        $workflow->add_operation($digital_expression_detection_operation);

        $digital_expression_detection_operation->lsf_queue($lsf_queue);
        $digital_expression_detection_operation->lsf_project($lsf_project);

        for my $key (keys %$params, qw/requestor sponsor user label output_dir/) {
            $workflow->connect_input(
                input_property => 'digital_expression_' . $key,
                destination => $digital_expression_detection_operation,
                destination_property => $key,
            );
        }

        $workflow->create_link(
            source => $alignment_operation,
            source_property => 'individual_alignment_results',
            destination => $digital_expression_detection_operation,
            destination_property => 'alignment_results',
        );

        $workflow->connect_output(
            source => $digital_expression_detection_operation,
            source_property => 'output_result',
            output_property => 'digital_expression_detection_result'
        );
    }


    # TopHat2 Specific Pipeline Steps
    if ($version_number >= 2) {
        if ($aligner_name eq 'tophat') {
            my $alignment_metrics_operation = Genome::WorkflowBuilder::Command->create(
                name => 'RnaSeq Alignment Metrics',
                command => 'Genome::Model::RnaSeq::Command::Tophat2AlignmentStats',
            );
            $workflow->add_operation($alignment_metrics_operation);
            $alignment_metrics_operation->lsf_queue($lsf_queue);
            $alignment_metrics_operation->lsf_project($lsf_project);

            $workflow->create_link(
                source => $alignment_operation,
                source_property => 'build_id',
                destination => $alignment_metrics_operation,
                destination_property => 'build_id'
            );
            $workflow->connect_output(
                source => $alignment_metrics_operation,
                source_property => 'result',
                output_property => 'alignment_stats_result'
            );
        }

        my $bam_qc_operation = Genome::WorkflowBuilder::Command->create(
            name => 'RnaSeq BamQc',
            command => 'Genome::Model::RnaSeq::Command::BamQc',
        );
        $workflow->add_operation($bam_qc_operation);
        $bam_qc_operation->lsf_queue($lsf_queue);
        $bam_qc_operation->lsf_project($lsf_project);

        $workflow->create_link(
            source => $alignment_operation,
            source_property => 'build_id',
            destination => $bam_qc_operation,
            destination_property => 'build_id'
        );
        $workflow->connect_output(
            source => $bam_qc_operation,
            source_property => 'result',
            output_property => 'bam_qc_result'
        );
        if ($run_splice_junction_summary) {
            my $splice_junction_operation = Genome::WorkflowBuilder::Command->create(
                name => 'RnaSeq Splice Junction Summary',
                command => 'Genome::Model::RnaSeq::Command::SpliceJunctionSummary',
            );
            $workflow->add_operation($splice_junction_operation);
            $splice_junction_operation->lsf_queue($lsf_queue);
            $splice_junction_operation->lsf_project($lsf_project);

            $workflow->create_link(
                source => $alignment_operation,
                source_property => 'build_id',
                destination => $splice_junction_operation,
                destination_property => 'build_id'
            );
            $workflow->connect_output(
                source => $splice_junction_operation,
                source_property => 'result',
                output_property => 'splice_junction_result'
            );
        }
    }

    # Picard
    my $picard_operation = Genome::WorkflowBuilder::Command->create(
        name => 'RnaSeq Picard Metrics',
        command => 'Genome::Model::RnaSeq::Command::PicardRnaSeqMetrics',
    );
    $workflow->add_operation($picard_operation);
    $picard_operation->lsf_queue($lsf_queue);
    $picard_operation->lsf_project($lsf_project);

    $workflow->create_link(
        source => $alignment_operation,
        source_property => 'build_id',
        destination => $picard_operation,
        destination_property => 'build_id'
    );

    # RefCov
    my $coverage_operation = Genome::WorkflowBuilder::Command->create(
        name => 'RnaSeq Coverage',
        command => 'Genome::Model::RnaSeq::Command::Coverage',
    );
    $workflow->add_operation($coverage_operation);
    $coverage_operation->lsf_queue($lsf_queue);
    $coverage_operation->lsf_project($lsf_project);

    $workflow->create_link(
        source => $alignment_operation,
        source_property => 'build_id',
        destination => $coverage_operation,
        destination_property => 'build_id'
    );

    # Cufflinks
    my $cufflinks_operation = Genome::WorkflowBuilder::Command->create(
        name => 'RnaSeq Cufflinks Expression',
        command => 'Genome::Model::RnaSeq::Command::Expression::Cufflinks',
    );
    $workflow->add_operation($cufflinks_operation);
    $cufflinks_operation->lsf_queue($lsf_queue);
    $cufflinks_operation->lsf_project($lsf_project);

    $workflow->create_link(
        source => $alignment_operation,
        source_property => 'build_id',
        destination => $cufflinks_operation,
        destination_property => 'build_id'
    );
    $workflow->connect_input(
        input_property => 'annotation_reference_transcripts_mode',
        destination => $cufflinks_operation,
        destination_property => 'annotation_reference_transcripts_mode'
    );

    $cufflinks_operation->parallel_by('annotation_reference_transcripts_mode');

    #Fusion Detection
    if($self->fusion_detector){
        my $operation_name = sprintf("RnaSeq Fusion Detection (%s %s)",
            $self->fusion_detector, $self->fusion_detector_version);
        my $fusion_detection_operation = Genome::WorkflowBuilder::Command->create(
            name => $operation_name,
            command => 'Genome::Model::RnaSeq::Command::DetectFusions',
        );
        $workflow->add_operation($fusion_detection_operation);
        $fusion_detection_operation->lsf_queue($lsf_queue);
        $fusion_detection_operation->lsf_project($lsf_project);

        $workflow->create_link(
            source => $alignment_operation,
            source_property => 'build_id',
            destination => $fusion_detection_operation,
            destination_property => 'build_id'
        );

        $self->connect_fusion_detector_to_input_connector($workflow, $fusion_detection_operation);

        #output connector
        $workflow->connect_output(
            source => $fusion_detection_operation,
            source_property => 'result',
            output_property => 'fusion_result'
        );

        ###Post Fusion Detection Annotation
        if($self->cancer_annotation_db){
            my $annotation_operation = Genome::WorkflowBuilder::Command->create(
                name => 'RnaSeq Fusion Detection Annotation',
                command => 'Genome::Model::RnaSeq::Command::AnnotateChimerascan',
            );
            $workflow->add_operation($annotation_operation);
            $annotation_operation->lsf_queue($lsf_queue);
            $annotation_operation->lsf_project($lsf_project);

            $workflow->create_link(
                source => $fusion_detection_operation,
                source_property => 'build_id',
                destination => $annotation_operation,
                destination_property => 'build_id',
            );

            $workflow->connect_input(
                input_property => 'cancer_annotation_db',
                destination => $annotation_operation,
                destination_property => 'cancer_annotation_db_id',
            );


            $workflow->connect_output(
                source => $annotation_operation,
                source_property => 'annotated_bedpe_file',
                output_property => 'annotated_bedpe_file',
            );
        }
    }

    # Define output connector results from coverage and expression
    $workflow->connect_output(
        source => $picard_operation,
        source_property => 'result',
        output_property => 'metrics_result'
    );
    $workflow->connect_output(
        source => $coverage_operation,
        source_property => 'result',
        output_property => 'coverage_result'
    );
    $workflow->connect_output(
        source => $cufflinks_operation,
        source_property => 'result',
        output_property => 'expression_result'
    );

    my $log_directory = $build->log_directory;
    $workflow->recursively_set_log_dir($log_directory);

    return $workflow;
}

sub connect_fusion_detector_to_input_connector {
    my $self = shift;
    my $workflow = shift;
    my $fusion_detection_operation = shift;

    my %input_properties = (
        'fusion_detector'           => 'detector_name',
        'fusion_detector_version'   => 'detector_version',
        'fusion_detector_params'    => 'detector_params',
    );

    while (my ($input, $destination) = each %input_properties) {
        $workflow->connect_input(
            input_property => $input,
            destination => $fusion_detection_operation,
            destination_property => $destination,
        );
    }

    return;
}

sub fusion_detection_inputs {
    my $self = shift;
    my $processing_profile = shift;

    return (
        fusion_detector => $processing_profile->fusion_detector,
        fusion_detector_version => $processing_profile->fusion_detector_version,
        fusion_detector_params => $processing_profile->fusion_detector_params,
    );
}

sub chimerascan_annotation_inputs {
    my $self = shift;
    my $build = shift;

    return (
        cancer_annotation_db => $build->cancer_annotation_db,
    );
}


sub params_for_alignment {
    my $self = shift;
    my @inputs = @_;

    my $reference_build = $self->reference_sequence_build;
    my $reference_build_id = $reference_build->id;

    my $read_aligner_params = $self->read_aligner_params || undef;

    if ($self->annotation_build) {
        my $annotation_build = $self->annotation_build;
        my $gtf_path = $annotation_build->annotation_file('gtf',$reference_build_id);
        unless (defined($gtf_path)) {
            die('There is no annotation GTF file defined for annotation_reference_transcripts build: '. $annotation_build->__display_name__);
        }

        # Test to see if this is version 1.4.0 or greater
        if (version->parse($self->read_aligner_version) >= version->parse('1.4.0')) {
            my $transcriptome_index_prefix = $annotation_build->annotation_file('',$reference_build_id);
            unless (-s $transcriptome_index_prefix .'.fa') {
                # TODO: We should probably lock until the first Tophat job completes creating the transriptome index
            }
            $read_aligner_params .= ' --transcriptome-index '. $transcriptome_index_prefix;
        }

        if ($read_aligner_params =~ /-G/) {
            die ('This processing_profile is requesting annotation_reference_transcripts \''. $annotation_build->__display_name__ .'\', but there seems to be a GTF file already defined in the read_aligner_params: '. $read_aligner_params);
        }
        if (defined($read_aligner_params)) {
            $read_aligner_params .= ' -G '. $gtf_path;
        } else {
            $read_aligner_params = ' -G '. $gtf_path;
        }
    }

    my %params = (
        instrument_data_id => [map($_->value_id, @inputs)],
        aligner_name => $self->read_aligner_name,
        reference_build_id => $reference_build_id || undef,
        aligner_version => $self->read_aligner_version || undef,
        aligner_params => $read_aligner_params,
        force_fragment => undef, #unused,
        trimmer_name => $self->read_trimmer_name || undef,
        trimmer_version => $self->read_trimmer_version || undef,
        trimmer_params => $self->read_trimmer_params || undef,
        picard_version => $self->picard_version || undef,
        samtools_version => $self->samtools_version || undef,
        filter_name => undef, #unused
        test_name => Genome::Config::get('software_result_test_name') || undef,
        bowtie_version => $self->bowtie_version
    );
    #$self->debug_message('The AlignmentResult parameters are: '. Data::Dumper::Dumper(%params));
    my @param_set = (\%params);
    return @param_set;
}

sub publication_description {
    my $self = shift;

    # TODO: use these, to dereive the values in the following two sections
    my $pp     = $self->processing_profile;
    my $refseq = $self->reference_sequence_build;
    my $annot  = $self->annotation_build;
    my @i      = $self->instrument_data;
    my $aligner_name = $self->read_aligner_name;

    # ensure we really only use one lane of data per library like we say we do
    my %libraries;
    for my $i (@i) {
        my $instdata_list = $libraries{$i->library_id} ||= [];
        if ($i->index_sequence) {
            die "the publication description is hard-coded to expect one lane of data per library";
        }
        push @$instdata_list, $i;
    }
    for my $library (keys %libraries) {
        my $i = $libraries{$library};
        if (@$i > 1) {
            die "the publication description is hard-coded to expect one lane of data per library";
        }
    }
    my $lane_count_summary = 'A single lane';

    # TODO: we must look this up from LIMS
    my $instrument = 'CHECKME HiSeq';
    my $chemistry  = 'CHECKME v3';
    my $lims_samtools_version = 'CHECKME 0.1.18';
    my $picard_version = 'CHECKME 1.4.6';

    # ensure that we are really on the build 37 reference
    my ($species, $alignment_ref);
    if ($refseq->id == 106942997) {
        $species = 'human';
        $alignment_ref = 'human reference genome (NCBI build 37)';
    }
    else {
        die "the publication description is hard-coded for human build 37 but got " . $refseq->id;
    }

    # ensure everything else we have hard-coded in the description still applies...
    my %expect = (
        read_aligner_name => $aligner_name,
        expression_name   => 'cufflinks',
    );
    for my $name (sort keys %expect) {
        my $expected_value = $expect{$name};
        my $actual_value = $self->$name;
        unless ($expected_value eq $actual_value) {
            die "publication description is hard-coded to expect that $name is '$expected_value', but got '$actual_value'";
        }
    }

    my $aligner_version   = $self->read_aligner_version;
    my $cufflinks_version = $self->expression_version;

    # TODO: update these to come from the model inputs and processing profile
    my $annotation_source = 'CHECKME the human Ensembl database (version 58) (REF)';
    my $bam_index_tool    = 'CHECKME samtools (v. 0.1.18)';
    my $bam_sort_tool     = 'CHECKME Picard (v.1.46)';

    my $file = __FILE__;
    my $line = __LINE__;

    my $desc = <<EOS;
RNA-seq analysis methods


$lane_count_summary of $instrument ($chemistry chemistry) was generated for each
Illumina RNA-seq library.  Reads were initially aligned to the $species
reference genome using Eland and stored as a BAM file.  These alignments
were used for basic quality assessment purposes only and no read filtering
was performed.  Mapping statistics for the BAM file were generated
using Samtools flagstat (v. $lims_samtools_version) (REF).
The BAM file was converted to FastQ using Picard (v.$picard_version) (REF)
and all reads were re-aligned to the $alignment_ref
using $aligner_name (v $aligner_version) (REF).  $aligner_name was run in default mode with
the following exceptions.  The --mate-inner-dist and --mate-std-dev
were estimated prior to run time using the Eland alignments described
above (elaborate) and specified at run time.  The '-G' option was used
to specify a GTF file for $aligner_name to generate an exon-exon junction
database to assist in the mapping of known junctions.  The transcripts
for this GTF were obtained from $annotation_source.  The
resulting $aligner_name BAM file was indexed by $bam_index_tool
and sorted by chromosome mapping position using $bam_sort_tool.
Transcript expression values were estimated by Cufflinks (v$cufflinks_version)
(REF) using default parameters with the following exceptions.  The Cufflinks
parameter '-G' was specified to force cufflinks to estimate expression
for known transcripts provided by the same GTF file that was supplied
to $aligner_name described above.  A second GTF containing only the
mitochondrial and ribosomal sequences was created and Cufflinks was
directed to ignore these regions using the '-M' mask option, to improve
overall robustness of transcript abundance estimates.  The variant
and corresponding gene expression status in the transcriptome were
determined for SNV positions identified as somatic in the WGS
tumor/normal data.  FPKM values were summarized to the gene level by
adding Cufflinks FPKMs from alternative transcripts of each Ensembl gene.
The variant allele frequencies were determined by counting reads
supporting reference and variant base counts using the Perl module
"Bio::DB::Sam".


Improve this description at line $line of file $file.

EOS

  $desc =~ s/\n(?!\n)/ /g;
  return $desc;
}

sub individual {
    my $self = shift;
    return $self->subject->individual;
}

1;

