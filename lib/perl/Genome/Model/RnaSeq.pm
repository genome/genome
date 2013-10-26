package Genome::Model::RnaSeq;

use strict;
use warnings;

use Genome;
use version;
use Genome::Utility::Text;

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
        },
        target_region_set_name => {
            is => 'Text',
            is_optional => 1,
            doc => 'limits the assignment of instrument data by default to only data with a matching TRSN'
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
        fusion_detection_strategy => {
            is_optional => 1,
            is => 'Text',
            doc => 'program, version and params to use for fusion detection ex: chimerascan 0.4.3 [-v]'
        },
        bowtie_version => {
            is_optional => 1,
            is => 'Text',
            doc => 'version of bowtie for tophat to use internally',
        }
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
    my $self = shift;
    my $build = shift;

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
            $params->{wrapper_build} = $build;
            $params->{wrapper_build_label} = $strategy_name . '_result';
            for my $key (keys %$params) {
                my $input_name = $strategy_name . '_' . $key;
                my $value = $params->{$key};
                push @inputs, $input_name => $value;
            }
        }
    }

    my %inputs = @inputs;

    return @inputs;
}

sub _resolve_workflow_for_build {
    # This is called by Genome::Model::Build::start()
    # Returns a Workflow::Operation
    # By default, builds this from stages(), but can be overridden for custom workflow.
    my $self = shift;
    my $build = shift;
    my $lsf_queue = shift; # TODO: the workflow shouldn't need this yet
    my $lsf_project = shift;

    if (!defined $lsf_queue || $lsf_queue eq '' || $lsf_queue eq 'inline') {
        $lsf_queue = 'apipe';
    }
    if (!defined $lsf_project || $lsf_project eq '') {
        $lsf_project = 'build' . $build->id;
    }

    my ($version_number) = $self->read_aligner_version =~ /^([\d+])/;

    my $processing_profile = $self->processing_profile;

    my $run_splice_junction_summary;

    if ($self->bedtools_version && $self->annotation_build) {
        my $annotation_build = $self->annotation_build;
        my ($ensembl_version) = split(/_/,$annotation_build->version);
        unless ($ensembl_version) {
            die('Failed to parse Ensembl version from annotation build: '. $annotation_build->id);
        }
        if ($ensembl_version >= 67) {
            $run_splice_junction_summary = 1;
        } else {
            $self->status_message('Skipping SpliceJunctionSummary for annotation build: '. $self->annotation_build);
        }
    }
    my $output_properties = ['coverage_result','expression_result','metrics_result'];
    push(@$output_properties, 'fusion_result') if $self->fusion_detection_strategy;
    push(@$output_properties, 'digital_expression_detection_result') if $self->digital_expression_detection_strategy;
    if ($version_number >= 2) {
        push(@$output_properties, 'bam_qc_result');
        push(@$output_properties, 'alignment_stats_result');
        if ($run_splice_junction_summary) {
            push(@$output_properties, 'splice_junction_result');
        }
    }

    my %inputs = $self->map_workflow_inputs($build);

    my $workflow = Workflow::Model->create(
        name => $build->workflow_name,
        input_properties => [keys %inputs],
        output_properties => $output_properties,
    );

    my $log_directory = $build->log_directory;
    $workflow->log_dir($log_directory);


    my $input_connector = $workflow->get_input_connector;
    my $output_connector = $workflow->get_output_connector;

    # Alignment
    my $alignment_operation = undef;
    if ($version_number < 2){
        $alignment_operation = $workflow->add_operation(
            name => 'RnaSeq Tophat Alignment',
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => 'Genome::Model::RnaSeq::Command::AlignReads::Tophat',
            )
        );
    }else{
        $alignment_operation = $workflow->add_operation(
            name => 'RnaSeq Alignment',
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => 'Genome::Model::RnaSeq::Command::AlignReads',
            )
        );
    }

    $alignment_operation->operation_type->lsf_queue($lsf_queue);
    $alignment_operation->operation_type->lsf_project($lsf_project);

    my $link = $workflow->add_link(
        left_operation => $input_connector,
        left_property => 'build_id',
        right_operation => $alignment_operation,
        right_property => 'build_id'
    );
   
    # Digital expression
    my $digital_expression_detection_operation = undef;
    if (my $strategy = $processing_profile->digital_expression_detection_strategy) {
        my ($class,$params) = $self->_parse_strategy($strategy, $build);

        $digital_expression_detection_operation = $workflow->add_operation(
            name => 'RnaSeq Digital Expression Detection',
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => $class . '::BuildStepWrapper',
            )
        );

        $digital_expression_detection_operation->operation_type->lsf_queue($lsf_queue);
        $digital_expression_detection_operation->operation_type->lsf_project($lsf_project);

        for my $key (keys %$params, qw/wrapper_build wrapper_build_label/) {
            my $link = $workflow->add_link(
                left_operation => $input_connector,
                left_property => 'digital_expression_' . $key,
                right_operation => $digital_expression_detection_operation,
                right_property => $key,
            );
        }
        
        $link = $workflow->add_link(
            left_operation => $alignment_operation,
            left_property => 'individual_alignment_results',
            right_operation => $digital_expression_detection_operation,
            right_property => 'alignment_results',
        );
        
        $link = $workflow->add_link(
            left_operation => $digital_expression_detection_operation,
            left_property => 'result',
            right_operation => $output_connector,
            right_property => 'digital_expression_detection_result'
        );
    }


    # TopHat2 Specific Pipeline Steps
    if ($version_number >= 2) {
        my $alignment_metrics_operation = $workflow->add_operation(
            name => 'RnaSeq Alignment Metrics',
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => 'Genome::Model::RnaSeq::Command::Tophat2AlignmentStats',
            )
        );
        $alignment_metrics_operation->operation_type->lsf_queue($lsf_queue);
        $alignment_metrics_operation->operation_type->lsf_project($lsf_project);

        $workflow->add_link(
            left_operation => $alignment_operation,
            left_property => 'build_id',
            right_operation => $alignment_metrics_operation,
            right_property => 'build_id'
        );
        $workflow->add_link(
            left_operation => $alignment_metrics_operation,
            left_property => 'result',
            right_operation => $output_connector,
            right_property => 'alignment_stats_result'
        );

        my $bam_qc_operation = $workflow->add_operation(
            name => 'RnaSeq BamQc',
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => 'Genome::Model::RnaSeq::Command::BamQc',
            )
        );
        $bam_qc_operation->operation_type->lsf_queue($lsf_queue);
        $bam_qc_operation->operation_type->lsf_project($lsf_project);
        
        $workflow->add_link(
            left_operation => $alignment_operation,
            left_property => 'build_id',
            right_operation => $bam_qc_operation,
            right_property => 'build_id'
        );
        $workflow->add_link(
            left_operation => $bam_qc_operation,
            left_property => 'result',
            right_operation => $output_connector,
            right_property => 'bam_qc_result'
        );
        if ($run_splice_junction_summary) {
            my $splice_junction_operation = $workflow->add_operation(
                name => 'RnaSeq Splice Junction Summary',
                operation_type => Workflow::OperationType::Command->create(
                    command_class_name => 'Genome::Model::RnaSeq::Command::SpliceJunctionSummary',
                )
            );
            $splice_junction_operation->operation_type->lsf_queue($lsf_queue);
            $splice_junction_operation->operation_type->lsf_project($lsf_project);
            
            $workflow->add_link(
                left_operation => $alignment_operation,
                left_property => 'build_id',
                right_operation => $splice_junction_operation,
                right_property => 'build_id'
            );
            $workflow->add_link(
                left_operation => $splice_junction_operation,
                left_property => 'result',
                right_operation => $output_connector,
                right_property => 'splice_junction_result'
            );
        }
    }
    
    # Picard
    my $picard_operation = $workflow->add_operation(
        name => 'RnaSeq Picard Metrics',
        operation_type => Workflow::OperationType::Command->create(
            command_class_name => 'Genome::Model::RnaSeq::Command::PicardRnaSeqMetrics',
        )
    );
    $picard_operation->operation_type->lsf_queue($lsf_queue);
    $picard_operation->operation_type->lsf_project($lsf_project);

    $workflow->add_link(
        left_operation => $alignment_operation,
        left_property => 'build_id',
        right_operation => $picard_operation,
        right_property => 'build_id'
    );

    # RefCov
    my $coverage_operation = $workflow->add_operation(
        name => 'RnaSeq Coverage',
        operation_type => Workflow::OperationType::Command->create(
            command_class_name => 'Genome::Model::RnaSeq::Command::Coverage',
        )
    );
    $coverage_operation->operation_type->lsf_queue($lsf_queue);
    $coverage_operation->operation_type->lsf_project($lsf_project);

    $workflow->add_link(
        left_operation => $alignment_operation,
        left_property => 'build_id',
        right_operation => $coverage_operation,
        right_property => 'build_id'
    );

    # Cufflinks
    my $cufflinks_operation = $workflow->add_operation(
        name => 'RnaSeq Cufflinks Expression',
        operation_type => Workflow::OperationType::Command->create(
            command_class_name => 'Genome::Model::RnaSeq::Command::Expression::Cufflinks',
        )
    );
    $cufflinks_operation->operation_type->lsf_queue($lsf_queue);
    $cufflinks_operation->operation_type->lsf_project($lsf_project);

    $workflow->add_link(
        left_operation => $alignment_operation,
        left_property => 'build_id',
        right_operation => $cufflinks_operation,
        right_property => 'build_id'
    );
    $workflow->add_link(
        left_operation => $input_connector,
        left_property => 'annotation_reference_transcripts_mode',
        right_operation => $cufflinks_operation,
        right_property => 'annotation_reference_transcripts_mode'
    );

    $cufflinks_operation->parallel_by('annotation_reference_transcripts_mode');

    #Fusion Detection
    if($self->fusion_detection_strategy){
        my ($detector, $version) = split(/\s+/, $self->fusion_detection_strategy);

        my $fusion_detection_operation = $workflow->add_operation(
            name => "RnaSeq Fusion Detection ($detector $version)",
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => 'Genome::Model::RnaSeq::Command::DetectFusions::' . Genome::Utility::Text::string_to_camel_case($detector,"-"),
            )
        );

        $fusion_detection_operation->operation_type->lsf_queue($lsf_queue);
        $fusion_detection_operation->operation_type->lsf_project($lsf_project);

        $workflow->add_link(
            left_operation => $alignment_operation,
            left_property => 'build_id',
            right_operation => $fusion_detection_operation,
            right_property => 'build_id'
        );

        #output connector
        $workflow->add_link(
            left_operation => $fusion_detection_operation,
            left_property => 'result',
            right_operation => $output_connector,
            right_property => 'fusion_result'
        );

    }

    # Define output connector results from coverage and expression
    $workflow->add_link(
        left_operation => $picard_operation,
        left_property => 'result',
        right_operation => $output_connector,
        right_property => 'metrics_result'
    );
    $workflow->add_link(
        left_operation => $coverage_operation,
        left_property => 'result',
        right_operation => $output_connector,
        right_property => 'coverage_result'
    );
    $workflow->add_link(
        left_operation => $cufflinks_operation,
        left_property => 'result',
        right_operation => $output_connector,
        right_property => 'expression_result'
    );

    return $workflow;
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
        aligner_name => 'tophat',
        reference_build_id => $reference_build_id || undef,
        aligner_version => $self->read_aligner_version || undef,
        aligner_params => $read_aligner_params,
        force_fragment => undef, #unused,
        trimmer_name => $self->read_trimmer_name || undef,
        trimmer_version => $self->read_trimmer_version || undef,
        trimmer_params => $self->read_trimmer_params || undef,
        picard_version => $self->picard_version || undef,
        samtools_version => undef, #unused
        filter_name => undef, #unused
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
        bowtie_version => $self->bowtie_version
    );
    #$self->status_message('The AlignmentResult parameters are: '. Data::Dumper::Dumper(%params));
    my @param_set = (\%params);
    return @param_set;
}

sub publication_description {
    my $self = shift;

    # TODO: use these, to dereive the values in the following two sections
    my $pp = $self->processing_profile;
    my $refseq = $self->reference_sequence_build;
    my $annot = $self->annotation_build;
    my @i = $self->instrument_data;

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
    my $chemistry = 'CHECKME v3';
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
        read_aligner_name => 'tophat',
        expression_name => 'cufflinks',
    );
    for my $name (sort keys %expect) {
        my $expected_value = $expect{$name};
        my $actual_value = $self->$name;
        unless ($expected_value eq $actual_value) {
            die "publication description is hard-coded to expect that $name is '$expected_value', but got '$actual_value'";
        }
    }

    my $tophat_version = $self->read_aligner_version;
    my $cufflinks_version = $self->expression_version;

    # TODO: update these to come from the model inputs and processing profile
    my $annotation_source = 'CHECKME the human Ensembl database (version 58) (REF)';
    my $bam_index_tool = 'CHECKME samtools (v. 0.1.18)';
    my $bam_sort_tool = 'CHECKME Picard (v.1.46)';

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
using Tophat (v $tophat_version) (REF).  Tophat was run in default mode with
the following exceptions.  The --mate-inner-dist and --mate-std-dev
were estimated prior to run time using the Eland alignments described
above (elaborate) and specified at run time.  The '-G' option was used
to specify a GTF file for Tophat to generate an exon-exon junction
database to assist in the mapping of known junctions.  The transcripts
for this GTF were obtained from $annotation_source.  The 
resulting tophat BAM file was indexed by $bam_index_tool
and sorted by chromosome mapping position using $bam_sort_tool. 
Transcript expression values were estimated by Cufflinks (v$cufflinks_version)
(REF) using default parameters with the following exceptions.  The Cufflinks 
parameter '-G' was specified to force cufflinks to estimate expression
for known transcripts provided by the same GTF file that was supplied
to TopHat described above.  A second GTF containing only the
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

1;

