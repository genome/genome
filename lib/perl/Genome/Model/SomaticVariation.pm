package Genome::Model::SomaticVariation;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticVariation {
    is  => 'Genome::ModelDeprecated',
    has_param => [
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
        tiering_version => {
            is => "Text",
            is_many => 0,
            valid_values => ['1','2', '3'],
            doc => "Version of tiering bed files to grab from the associated annotation build",
        },
        loh_version => {
            is => "Text",
            is_many => 0,
            is_optional => 1,
            doc => "Version of loh used for this build. If none is specified, no loh detection will be done.",
        },
        transcript_variant_annotator_version => {
            doc => 'Version of the "annotate transcript-variants" tool to run during the annotation step',
            is_optional => 1,
            default_value => Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version,
            valid_values => [ 0,1,2,3,4 ],
        },
        get_regulome_db => {
            doc => "Get the regulome-db annotation for the snvs list",
            is_optional => 1,
            default_value => 0,
            is => "Boolean",
        },
        filter_previously_discovered_variants => {
            doc => 'Should variants be divided into previously discovered and novel variants',
            default_value => 0,
            is => 'Boolean',
        },
        vcf_annotate_dbsnp_info_field_string => {
            doc => 'String that indicates which dbsnp info fields should be included in the annotated vcf',
            is => 'Text',
            default_value => "NO_INFO",
        },
        required_snv_callers => {
            is => 'Number',
            doc => "Number of independent algorithms that must call a SNV in order to be included in the final report. If set to 1 (default), all calls are used",
            is_optional => 1,
            default_value => 1,
        },
        tiers_to_review => {
            is => 'String',
            doc => "comma-separated list of tiers to include in the final report. (e.g. 1,2,3 will create bed files with tier1, tier2, and tier3)",
            is_optional => 1,
            default_value => 1,
        },
   ],
    has => [
        tumor_model_id => {
            via => 'tumor_model',
            to => 'id',
        },
        normal_model_id => {
            via => 'normal_model',
            to => 'id',
        },
        annotation_build_id => {
            via => 'annotation_build',
            to => 'id',
        },
        previously_discovered_variations_build_id => {
            via => 'previously_discovered_variations',
            to => 'id',
            doc => 'previous variants genome feature set to screen somatic mutations against',
        },
        previously_discovered_variations_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            via => '__self__',
            to => 'previously_discovered_variations',
        },
        force => {
            is => 'Boolean',
            is_optional => 1,
            is_many => 0,
            default => 0,
            doc => 'Allow creation of somatic variation models where --tumor_model and --normal_model do not have matching Genome::Individuals',
        },
    ],
    has_input => [
        tumor_model => {
            is => 'Genome::Model',
        },
        normal_model => {
            is => 'Genome::Model',
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
        },
        previously_discovered_variations => {
            is => 'Genome::Model::Build::ImportedVariationList',
        },
        regulatory_annotations => {
            is => 'Genome::FeatureList',
            is_many => 1,
            is_optional => 1,
        },
    ],
    has_optional => [
        experimental_subject => {
            is => 'Genome::Sample',
            via => 'tumor_model',
            to => 'subject'
        },
        control_subject => {
            is => 'Genome::Sample',
            via => 'normal_model',
            to => 'subject'
        },
    ],
};

sub help_synopsis_for_create_profile {
    my $self = shift;
    return <<"EOS"
  Complete Examples:

    genome processing-profile create somatic-variation \
      --name 'unfiltered sniper with breakdancer' \
      --snv-detection-strategy 'sniper 0.7.3 [ -q 1 -Q 15 ] intersect samtools r599' \
      --indel-detection-strategy   '(sniper 0.7.3 [-q 1 -Q 15] filtered by library-support v1) union (samtools r599  intersect pindel 0.1)' \
      --sv-detection-strategy 'breakdancer 2010_06_24  filtered by tigra-assembly v1'

    genome processing-profile create somatic-variation \
      --name 'filtered sniper with breakdancer' \
      --snv-detection-strategy '(sniper 0.7.3 [-q 1 -Q 15] filtered by loh v1, somatic-score-mapping-quality v1 [-min_somatic_quality 40 -min_mapping_quality 40]) intersect samtools r599'  \
      --indel-detection-strategy 'sniper 0.7.3 [-q 1 -Q 15] filtered by library-support v1' \
      --sv-detection-strategy 'breakdancer 2010_06_24 filtered by tigra-assembly v1'

  Example Strategies usable for SNVs, indels, SVs, or combinations:

    'sniper 0.7.3 [-q 1 -Q 15]'
    # Detect with sniper version 0.7.3 with the parameters "-q 1 -Q 15".
    # works for SNVs, indels   

    'sniper 0.7.3 [-q 1 -Q 15] filtered by loh v1 '
    # Detect snvs or indels with sniper version 0.7.3 with the listed parameters and filter the results by running the "loh" filter version "v1".

    'sniper 0.7.3 [-q 1 -Q 15] filtered by loh v1, somatic-score-mapping-quality v1 [-min_somatic_quality 40:-min_mapping_quality 40] intersect samtools r599'  
    # Detect snvs and/or indels with the above as follows:
    # 1) Run sniper version 0.7.3 with parameters
    # 2) Filter the results by running the loh filter version v1
    # 3) Further filter results by running the somatic-score-mapping-quality filter version v1 with parameters.
    # 4) Run samtools version r599 (or steal previous results) 
    # 5) Intersect 3 & 4 

    'sniper 0.7.3 [-q 1 -Q 15] union (samtools r599  intersect pindel v1)'
    # Detect indels with: 
    # 1) Run sniper version 0.7.3 with the listed parameters. 
    # 2) Run samtools version r599 
    # 3) Run pindel version v1
    # 4) Intersect 2 and 3
    # 5) Union 1 and 4.
EOS
}

sub help_detail_for_create_profile {
    return <<EOS
  For a detailed explanation of how to writing a variant detection strategy, see:
    perldoc Genome::Model::Tools::DetectVariants2::Strategy;
EOS
}

# XXX This looks useful... but is it ever delivered to the user?
sub help_manual_for_create_profile {
    return <<EOS

  EXAMPLES

    Strategies usable for SNVs, indels, SVs, or combinations:

    'sniper 0.7.3 [-q 1 -Q 15]'
    # Detect with sniper version 0.7.3 with the parameters "-q 1 -Q 15".
    # works for SNVs, indels   

    'sniper 0.7.3 [ -q 1 -Q 15 ] filtered by loh v1 '
    # Detect snvs or indels with sniper version 0.7.3 with the listed parameters and filter the results by running the "loh" filter version "v1".

    'sniper 0.7.3 [ -q 1 -Q 15 ] filtered by loh v1 , somatic-score-mapping-quality v1 [-min_somatic_quality 40:-min_mapping_quality 40] intersect samtools r599'  
    # Detect snvs and/or indels with the above as follows:
    # 1) Run sniper version 0.7.3 with parameters
    # 2) Filter the results by running the loh filter version v1, i
    # 3) Further filter results and then the somatic-score-mapping-quality filter version v1 with parameters.
    # 4) Run samtools version r599 (or steal previous results) 
    # 5) Intersect 3 & 4 

    'sniper 0.7.3 [ -q 1 -Q 15 ] union (samtools r599  intersect pindel v1 )'
    # Detect indels with: 
    # 1) Run sniper version 0.7.3 with the listed parameters. 
    # 2) Run samtools version r599 
    # 3) Run pindel version v1
    # 4) Intersect 2 and 3
    # 5) Union 1 and 4.

    'sniper 0.7.3 [ -q 1 -Q 15 ]' 
    # Detect snvs or indels or both with sniper version 0.7.3 with the listed parameters. 
    # This expression can be set as an snv detection strategy or an indel detection strategy, 
    # and if both are set to the same value sniper will run just once to do both.

    'breakdancer 2010_06_24 ' 
    # Detect structural variation with breakdancer version 2010_06_24.

    'sniper 0.7.3 [ -q 1 -Q 15 ] intersect samtools r599 '
    # Detect snvs: Intersect the results of sniper version 0.7.3 with parameters and samtools version r599.

    'sniper 0.7.3 [ -q 1 -Q 15 ] filtered by library-support v1 '
    # Detect indels using sniper version 0.7.3 with parameters and filter the results with the library-support filter version v1

    'breakdancer 2010_06_24  filtered by tigra-assembly v1 '
    # Detect structural variations using breakdancer version 2010_06_24 and filter the results by applying the tigra-assembly filter version v1

    'sniper 0.7.3 [ -q 1 -Q 15 ] filtered by library-support v1 '
    # Detect indels using sniper version 0.7.3 with parameters and filter the results with the library-support filter version v1

    'breakdancer 2010_06_24  filtered by tigra-assembly v1 '
    # Detect structural variations using breakdancer version 2010_06_24 and filter the results by applying the tigra-assembly filter version v1

  EXPLANATION

    A strategy consists of the following:
    detector-name version [ params ] filtered by filter-name version [ params ],filter-name version [ params ] ...

    * Detector-name is the name of the variant detector as it follows "gmt detect-variants2". For example, "sniper" would reference the tool located at "gmt detect-variants2 sniper".

    * In the same way, filter-name is the name of the filter as it follows "gmt detect-variants2 filter". For example, "loh" would reference the tool located at gmt detect-variants2 filter loh".

    * Version is a version number that pertains to that detector or filter specifically. For sniper this might be "0.7.3". For samtools this might be "r599".
        Many filters are not currently versioned, but may be in the future. In these cases "v1" should be used to denote version 1.

    * The parameter list is a list of all parameters to be passed to the detector or filter and will be specific to that tool. It is passed as a single string and is optional.

    * Filtered by may contain any number of complete filter specifications (separated by commas), including 0. Each filter must be a complete list of name, version, and an optional param list.

    --- Unions and intersections ---

    * Variant detectors can be intersected or unioned with each other to create variant lists which utilize more than one variant detector. In either case, all variant detectors will be run individually and then processed together.
    * An intersection will run both detectors and then produce a final list of variants that represents where both detectors agree on both the position and the call.
    * A union will run both detectors and then produce a final list of variants that represents every call that both the detectors made, regardless of agreement.
    * Parenthesis may also be utilized around pairs of detectors to specify logical order of operation.

    --- Examples of union and intersection ---
    --snv-detection-strategy 'sniper 0.7.3 [-q 1 -Q 15] intersect samtools r599
    This represents the desire to run version 0.7.3 of sniper with the above parameter list and version r599 of samtools with no parameters and intersect the results.
    Both detectors will be run and the final variant list will represent all variants which were called by both detectors.

    --snv-detection-strategy 'sniper 0.7.3 [-q 1 -Q 15] union (samtools r599 intersect pindel v1)
    This represents the desire to run version 0.7.3 of sniper with the above parameters, version r599 of samtools with no parameters, and version v1 of pindel with no parameters.
    Due to the parenthesis, the results of pindel and samtools will first be intersected and then that result will be unioned with the variant calls from sniper.
    In plain language, the resulting set will be any variants that either a) sniper called or b) pindel and samtools both called and agreed on.

EOS
}

sub create {
    my $class  = shift;
    my $bx = $class->define_boolexpr(@_);
    my %params = $bx->params_list;

    my $dbsnp_flag = 1 if $params{previously_discovered_variations_build} or $params{previously_discovered_variations_build_id};

    my $tumor_model  = $params{tumor_model} || Genome::Model->get($params{tumor_model_id});
    my $normal_model =  $params{normal_model}  || Genome::Model->get($params{normal_model_id});;
    my $annotation_build = $params{annotation_build} || Genome::Model::Build->get($params{annotation_build_id});
    my $previously_discovered_variations_build = $params{previously_discovered_variations_build} || Genome::Model::Build->get($params{previously_discovered_variations_build_id})
        if $dbsnp_flag;

    unless($tumor_model) {
        $class->error_message('No tumor model provided.');
        return;
    }

    unless($normal_model) {
        $class->error_message('No normal model provided.');
        return;
    }

    unless($annotation_build) {
        $class->error_message('No annotation build provided.');
        return;
    }

    if ($dbsnp_flag) {
        unless($previously_discovered_variations_build) {
            $class->error_message('No previous variants build provided.');
            return;
        }
    }

    my $tumor_subject  = $tumor_model->subject;
    my $normal_subject = $normal_model->subject;

    if ($tumor_subject->can('source') and $normal_subject->can('source')) {

        my $tumor_source  = $tumor_subject->source;
        my $normal_source = $normal_subject->source;

        unless ($tumor_source eq $normal_source) {
            my $tumor_common_name  = $tumor_source->common_name  || "unknown";
            my $normal_common_name = $normal_source->common_name || "unknown";
            my $message = "Tumor model and normal model samples do not come from the same individual.  Tumor common name is $tumor_common_name. Normal common name is $normal_common_name.";
            if (defined $params{force} and $params{force} == 1) {
                $class->warning_message($message);
            }
            else{
                die $class->error_message($message . " Use --force to allow this anyway.");
            }
        }
        $params{subject_id} = $tumor_subject->id;
        $params{subject_class_name} = $tumor_subject->class;

    } 
    else {
        $class->error_message('Unexpected subject for tumor or normal model!');
        return;
    }

    my $self = $class->SUPER::create(%params);

    unless ($self){
        $class->error_message('Error in model creation');
        return;
    }

    unless($self->tumor_model) {
        $self->error_message('No tumor model on model!' );
        $self->delete;
        return;
    }

    unless($self->normal_model) {
        $self->error_message('No normal model on model!');
        $self->delete;
        return;
    }

    unless($self->annotation_build) {
        $self->error_message('No annotation build on model!' );
        $self->delete;
        return;
    }

    if ($dbsnp_flag) {
        unless($self->previously_discovered_variations_build) {
            $self->error_message('No previous variants build on model!');
            $self->delete;
            return;
        }
    }

    return $self;
}

sub _input_differences_are_ok {
    my $self = shift;
    my @inputs_not_found = @{shift()};
    my @build_inputs_not_found = @{shift()};

    return unless scalar(@inputs_not_found) == 2 and scalar(@build_inputs_not_found) == 2;

    my $input_sorter = sub { $a->name cmp $b->name };

    @inputs_not_found = sort $input_sorter @inputs_not_found;
    @build_inputs_not_found = sort $input_sorter @build_inputs_not_found;

    #these are expected to differ and no new build is needed as long as the build pointed to is the latest for the model
    for(0..$#inputs_not_found) {
        return unless $inputs_not_found[$_]->value && $inputs_not_found[$_]->value->isa('Genome::Model');
        return unless $inputs_not_found[$_]->value->last_complete_build eq $build_inputs_not_found[$_]->value;
    }

    return 1;
}

sub _initialize_build {
    my($self,$build) = @_;
    return 1;
}

sub _resolve_workflow_for_build {
    my $self = shift;
    my $build = shift;

    my $operation = Workflow::Operation->create_from_xml(__FILE__ . '.xml');

    my $log_directory = $build->log_directory;
    $operation->log_dir($log_directory);

    #I think this ideally should be handled 
    $operation->name($build->workflow_name);

    return $operation;
}

sub map_workflow_inputs {
    my $self = shift;
    my $build = shift;

    my @inputs = ();

    # Verify the somatic model
    my $model = $build->model;

    unless ($model) {
        $self->error_message("Failed to get a model for this build!");
        die $self->error_message;
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

    my $tumor_bam;
    if ($tumor_build->isa('Genome::Model::Build::RnaSeq')) {
        $tumor_bam = $tumor_build->merged_alignment_result->bam_file;
    } else {
        $tumor_bam = $tumor_build->whole_rmdup_bam_file;
    }
    unless (-e $tumor_bam) {
        $self->error_message("Tumor bam file $tumor_bam does not exist!");
        die $self->error_message;
    }
    my $normal_bam;
    if ($normal_build->isa('Genome::Model::Build::RnaSeq')) {
        $normal_bam = $normal_build->merged_alignment_result->bam_file;
    } else {
        $normal_bam = $normal_build->whole_rmdup_bam_file;
    }
    unless (-e $normal_bam) {
        $self->error_message("Normal bam file $normal_bam does not exist!");
        die $self->error_message;
    }

    push @inputs, build_id => $build->id;
    push @inputs, annotator_version => $self->transcript_variant_annotator_version;
    push @inputs, regulatory_annotations => [$self->regulatory_annotations];
    push @inputs, get_regulome_db => $self->get_regulome_db;
    push @inputs, required_snv_callers => $self->required_snv_callers;
    push @inputs, tiers_to_review => $self->tiers_to_review;

    return @inputs;
}

sub __profile_errors__ {
    my ($self, $pp) = @_;

    my @errors;
    if ($pp->snv_detection_strategy) {
        my $snv_strat = Genome::Model::Tools::DetectVariants2::Strategy->get($pp->snv_detection_strategy);
        push @errors, $snv_strat->__errors__;
    }
    if ($pp->sv_detection_strategy) {
        my $sv_strat = Genome::Model::Tools::DetectVariants2::Strategy->get($pp->sv_detection_strategy);
        push @errors, $sv_strat->__errors__;
    }
    if ($pp->indel_detection_strategy) {
        my $indel_strat = Genome::Model::Tools::DetectVariants2::Strategy->get($pp->indel_detection_strategy);
        push @errors, $indel_strat->__errors__;
    }
    return @errors
}

1;
