package Genome::Model::SomaticCapture;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticCapture {
    is  => 'Genome::ModelDeprecated',
    has_param => [
        only_tier_1 => {
            doc => "If set to true, the pipeline will skip ucsc annotation and produce only tier 1 snps",
        },
        min_mapping_quality => {
            doc => "minimum average mapping quality threshold for high confidence call",
        },
        min_somatic_quality => {
            doc => "minimum somatic quality threshold for high confidence call",
        },
        skip_sv => {
            doc => "If set to true, the pipeline will skip structural variation detection",
        },
        transcript_variant_annotator_version => {
            doc => 'Version of the "annotate transcript-variants" tool to run during the annotation step',
            is_optional => 1,
            default_value => Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version,
            valid_values => [ 0,1,2,3],#Genome::Model::Tools::Annotate::TranscriptVariants->available_versions ],
        },
    ],
    has_optional => [
        tumor_model_links => {
            is => 'Genome::Model::Link',
            reverse_as => 'to_model',
            where => [ role => 'tumor'],
            is_many => 1,
        },
        tumor_model => {
            is => 'Genome::Model',
            via => 'tumor_model_links',
            to => 'from_model',
        },
        tumor_model_id => {
            is => 'Integer',
            via => 'tumor_model',
            to => 'id',
        },
        normal_model_links => {
            is => 'Genome::Model::Link',
            reverse_as => 'to_model',
            where => [ role => 'normal'],
            is_many => 1,
        },
        normal_model => {
            is => 'Genome::Model',
            via => 'normal_model_links',
            to => 'from_model',
        },
        normal_model_id => {
            is => 'Integer',
            via => 'normal_model',
            to => 'id',
        },
    ],
};

sub create {
    my $class = shift;
    my %params = @_;

    my $tumor_model_id = delete $params{tumor_model_id};
    my $tumor_model = delete $params{tumor_model};
    my $normal_model_id = delete $params{normal_model_id};
    my $normal_model = delete $params{normal_model};

    unless($tumor_model) {
        $tumor_model = Genome::Model->get($tumor_model_id);

        unless($tumor_model) {
            $class->error_message('Could not find tumor model.' );
            return;
        }
    }

    unless($normal_model) {
        $normal_model = Genome::Model->get($normal_model_id);

        unless($normal_model) {
            $class->error_message('Could not find normal model.');
            return;
        }
    }

    my $tumor_subject = $tumor_model->subject;
    my $normal_subject = $normal_model->subject;

    if(($tumor_subject->can('source') || $tumor_subject->can('sample')) and ($normal_subject->can('source') || $normal_subject->can('sample'))) {

        my $tumor_source;
        if($tumor_subject->can('source')) {
            $tumor_source = $tumor_subject->source;
        } else {
            $tumor_source = $tumor_subject->sample->source;
        }

        my $normal_source;
        if($normal_subject->can('source')) {
            $normal_source = $normal_subject->source;
        } else {
            $normal_source = $normal_subject->sample->source;
        }

        if($tumor_source eq $normal_source) {
            my $subject = $tumor_source;

            #Set up other parameters for call to parent execute()
            $params{subject_id} = $subject->id;
            $params{subject_class_name} = $subject->class;
        } else {
            $class->error_message('Tumor and normal samples are not from same source!');
            return;
        }
    } else {
        $class->error_message('Unexpected subject for tumor or normal model!');
        return;
    }

    my $self = $class->SUPER::create(%params);

    $self->add_from_model(from_model => $normal_model, role => 'normal');
    $self->add_from_model(from_model => $tumor_model, role => 'tumor');

    return $self;
}

sub get_all_objects {
    my $self = shift;

    my @objects = $self->SUPER::get_all_objects(@_);
    my @validations = Genome::Model::VariantValidation->get(model_id=>$self->id);
    push @objects, @validations;
    return @objects;
}

sub _resolve_workflow_for_build {
    my $self = shift;
    my $build = shift;

    my $operation = Workflow::Operation->create_from_xml(__FILE__ . '.xml');

    my $log_directory = $build->log_directory;
    $operation->log_dir($log_directory);
    $operation->name($build->workflow_name);

    return $operation;
}

sub _initialize_build {
    my $self = shift;
    my $build = shift;

    my $tumor_model = $self->tumor_model;
    unless($tumor_model) {
        $self->error_message('No tumor model found.');
        return;
    }
    my $normal_model = $self->normal_model;
    unless($normal_model) {
        $self->error_message('No normal model found.');
        return;
    }

    my $tumor_build = $tumor_model->last_succeeded_build;
    unless($tumor_build) {
        $self->error_message('No succeeded build found for the tumor model.');
        return;
    }
    my $normal_build = $normal_model->last_succeeded_build;
    unless($normal_build) {
        $self->error_message('No succeeded build found for the normal model.');
        return;
    }

    $build->add_from_build(from_build => $tumor_build, role => 'tumor');
    $build->add_from_build(from_build => $normal_build, role => 'normal');

    # Check that the annotator version param is sane before doing the build
    my $annotator_version;
    my $worked = eval {
        $annotator_version = $self->transcript_variant_annotator_version;
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

sub map_workflow_inputs {
    my $self = shift;
    my $build = shift;

    my @inputs = ();

    push @inputs,
        build_id => $build->id;

    if ($build->data_directory) {
        unless (-d $build->data_directory) {
            $self->error_message("Data directory " . $build->data_directory . " does not exist. Please create it.");
            return 0;
        }

        push @inputs,
            data_directory => $build->data_directory;

        my %default_filenames = $self->default_filenames;
        for my $param (keys %default_filenames) {
            # set a default param if one has not been specified
            my $default_filename = $default_filenames{$param};
            push @inputs, $param => join('/', $build->data_directory, $default_filename );
        }
    }

    my $tumor_build = $build->tumor_build;
    my $normal_build = $build->normal_build;

    unless ($tumor_build) {
        $self->error_message("Failed to get a tumor_build associated with this somatic capture build!");
        die $self->error_message;
    }

    unless ($normal_build) {
        $self->error_message("Failed to get a normal_build associated with this somatic capture build!");
        die $self->error_message;
    }

    my $data_directory = $build->data_directory;
    unless ($data_directory) {
        $self->error_message("Failed to get a data_directory for this build!");
        die $self->error_message;
    }

    my $tumor_bam = $tumor_build->whole_rmdup_bam_file;
    unless (-e $tumor_bam) {
        $self->error_message("Tumor bam file $tumor_bam does not exist!");
        die $self->error_message;
    }

    my $normal_bam = $normal_build->whole_rmdup_bam_file;
    unless (-e $normal_bam) {
        $self->error_message("Normal bam file $normal_bam does not exist!");
        die $self->error_message;
    }

    # Get the snp file from the tumor and normal models
    my $tumor_snp_file = $tumor_build->snv_file;
    unless (-e $tumor_snp_file) {
        $self->error_message("Tumor snp file $tumor_snp_file does not exist!");
        die $self->error_message;
    }
    my $normal_snp_file = $normal_build->snv_file;
    unless (-e $normal_snp_file) {
        $self->error_message("Normal snp file $normal_snp_file does not exist!");
        die $self->error_message;
    }

    ## Set a default reference transcript annotator version ##

    my $annotation_reference_transcripts = "NCBI-human.combined-annotation/54_36p_v2";
    my $ref_id = $tumor_build->model->reference_sequence_build->id;

    ## Human build 37 ##    
    if($ref_id == 106942997 || $ref_id == 102671028)
    {
        $annotation_reference_transcripts = "NCBI-human.combined-annotation/57_37b";
    }
    elsif($ref_id == 104420993 || $ref_id == 103785428)
    {
        $annotation_reference_transcripts = "NCBI-mouse.combined-annotation/54_37g_v2";
    }

    push @inputs,
        tumor_bam_file => $tumor_bam,
        normal_bam_file => $normal_bam,
        tumor_snp_file => $tumor_snp_file,
        normal_snp_file => $normal_snp_file; 

    my $reference_fasta = $tumor_build->model->reference_sequence_build->full_consensus_path('fa');

    # Set (hardcoded) defaults for tools that have defaults that do not agree with somatic pipeline
    push @inputs,
        skip_if_output_present => 1,
        lookup_variants_report_mode => 'novel-only',
        lookup_variants_filter_out_submitters => "SNP500CANCER,OMIMSNP,CANCER-GENOME,CGAP-GAI,LCEISEN,ICRCG",
        annotate_no_headers => 1,
        transcript_annotation_filter => 'top',
        only_tier_1_indel => 1,
        normal_indelpe_data_directory => join('/', $build->data_directory, "normal_indelpe_data" ),
        tumor_indelpe_data_directory => join('/', $build->data_directory, "tumor_indelpe_data" ),
        reference_fasta => $reference_fasta,
        annotation_reference_transcripts => $annotation_reference_transcripts,
        prepend_chr => 0;

    # Set values from processing profile parameters
    push @inputs,
        only_tier_1 => (defined($self->only_tier_1)? $self->only_tier_1 : 0),
        skip_sv => (defined($self->skip_sv)? $self->skip_sv : 0),
        transcript_variant_annotator_version => (defined($self->transcript_variant_annotator_version)? $self->transcript_variant_annotator_version : ':'),
        min_mapping_quality => (defined($self->min_mapping_quality) ? $self->min_mapping_quality : 40),
        min_somatic_quality => (defined($self->min_somatic_quality) ? $self->min_somatic_quality : 40);

    return @inputs;
}

sub default_filenames{
    my $self = shift;

    my %default_filenames = (
        breakdancer_working_directory       => 'breakdancer/',
        sniper_working_directory            => 'sniper/',


        ## glfSomatic Output Files ##
        sniper_snp_output_adaptor           => 'somaticSniper.output.snp.adaptor',
        sniper_snp_output_filter            => 'somaticSniper.output.snp.filter',
        sniper_snp_output_filter_hc         => 'somaticSniper.output.snp.filter.hc',
        sniper_snp_output_filter_hc_somatic => 'somaticSniper.output.snp.filter.hc.somatic',
        sniper_snp_output_filter_hc_loh     => 'somaticSniper.output.snp.filter.hc.loh',

        ## Files formatted for annotation ##

        adaptor_output_indel                => 'somaticSniper.output.indel.formatted',
        filter_indel_output                 => 'somaticSniper.output.indel.formatted.filter',

        ## Varscan Output Files ##
        varscan_snp_output                  => 'varScan.output.snp',
#        varscan_snp_output_filter           => 'varScan.output.snp.filter',
        varscan_indel_output                => 'varScan.output.indel',

        varscan_adaptor_snp                 => 'varScan.output.snp.formatted',
        varscan_adaptor_indel               => 'varScan.output.indel.formatted',
        varscan_snp_germline                => 'varScan.output.snp.formatted.Germline',
        varscan_snp_loh                     => 'varScan.output.snp.formatted.LOH',
        varscan_snp_somatic                 => 'varScan.output.snp.formatted.Somatic.hc',

        varscan_indel_germline              => 'varScan.output.indel.formatted.Germline',
        varscan_indel_loh                   => 'varScan.output.indel.formatted.LOH',
        varscan_indel_somatic               => 'varScan.output.indel.formatted.Somatic',

        ## GATK Files ##

        gatk_output                         => 'gatk.output.indel',
        gatk_output_formatted               => 'gatk.output.indel.formatted',
        gatk_output_somatic                 => 'gatk.output.indel.formatted.Somatic',

        ## Combined glfSomatic+Varscan Output files ##
        merged_snp_output                   => 'merged.somatic.snp',            ## Generated from merge-variants of somaticSniper and varScan
        merged_indel_output                 => 'merged.somatic.indel',          ## Generated from merge-variants of somaticSniper and varScan ##
        merged_indel_output_filter          => 'merged.somatic.indel.filter',          ## The fp-filtered list of indels ##
        merged_indel_output_filter_rem          => 'merged.somatic.indel.filter.removed',          ## The fp-filter removed list of indels ##

        ## Strand Filtering and Lookup Variants Files ##
        merged_snp_output_filter        => 'merged.somatic.snp.filter',
        merged_snp_output_filter_fail   => 'merged.somatic.snp.filter.removed',
        merged_snp_output_novel         => 'merged.somatic.snp.filter.novel',


        ## Annotation output files ##
        annotate_output_snp                 => 'annotation.somatic.snp.transcript',
        ucsc_output                         => 'annotation.somatic.snp.ucsc',
        ucsc_unannotated_output             => 'annotation.somatic.snp.unannot-ucsc',

        annotate_output_indel                 => 'annotation.somatic.indel.transcript',

        ## Tiered SNP and indel files (all confidence) ##

        tier_1_snp_file                     => 'merged.somatic.snp.filter.novel.tier1',
        tier_2_snp_file                     => 'merged.somatic.snp.filter.novel.tier2',
        tier_3_snp_file                     => 'merged.somatic.snp.filter.novel.tier3',
        tier_4_snp_file                     => 'merged.somatic.snp.filter.novel.tier4',
        tier_1_indel_file                   => 'merged.somatic.indel.filter.tier1',

        ## Tiered SNP/indel files (high and highest conf ) ##

        tier_1_snp_file_high                => 'merged.somatic.snp.filter.novel.tier1.hc',
        tier_1_snp_file_highest             => 'merged.somatic.snp.filter.novel.tier1.gc',
        tier_1_indel_file_high              => 'merged.somatic.indel.filter.tier1.hc',
        tier_1_indel_file_highest           => 'merged.somatic.indel.filter.tier1.gc',

        ## Breakdancer and Copy Number files ##

        copy_number_output                  => 'copy_number.csv',
        circos_graph                        => 'circos_graph',
        variant_report_output               => 'cancer_report.html', 
        file_summary_report_output          => 'file_summary_report.html',


        upload_variants_snp_1_output        => 'upload-variants.snp_1.out',
        upload_variants_snp_2_output        => 'upload-variants.snp_2.out',
        upload_variants_indel_output        => 'upload-variants.indel.out',
    );

    return %default_filenames;
}

1;
