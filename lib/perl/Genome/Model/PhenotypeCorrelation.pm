package Genome::Model::PhenotypeCorrelation;

use strict;
use warnings;
use Genome;
use List::Util qw(reduce);
use Genome::File::Vcf::Reader;
use Math::Complex;
use File::chdir;
use File::Basename qw( fileparse );
use Digest::MD5;
use Carp qw(confess);

class Genome::Model::PhenotypeCorrelation {
    is => 'Genome::ModelDeprecated',
    doc => "genotype-phenotype correlation of a population group",
    has_param => [
        alignment_strategy => {
            is => "Text",
            is_many => 0,
            doc => "Strategy align sequence reads.",
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
        roi_wingspan => {
            is => 'Text',
            doc => 'Area to include before and after ROI regions',
            is_optional => 1,
        },
        group_samples_for_genotyping_by => {
            is => "Text",
            is_many => 0,
            is_optional => 1,
            #default_value => 'each',
            valid_values => ['each', 'trio', 'all'],
            doc => "group samples together when genotyping, using this attribute, instead of examining genomes independently (use \"all\" or \"trio\")",
        },
        trait_type => {
            is => "Text",
            is_many => 0,
            is_optional => 0,
            valid_values => [qw( quantitative case-control mendelian none)],
            doc => "Type of trait under analysis",
        },
        cohort_type => {
            is => "Text",
            is_many => 0,
            is_optional => 0,
            valid_values => [qw( unrelated family-based )],
            doc => "Type of cohort under study",
        },
        joinx_version => {
            is => "Text",
            doc => "Version of joinx to use",
            is_optional => 0,
            default => 1.6,
        },
        maximum_maf => {
            is => "Number",
            doc => "Maximum minor allele frequency to include in any downstream burden tests",
            is_optional => 1,
            default => 0.01,
        },
    ],
    has_input => [
        reference_sequence_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'the reference sequence against which alignment and variant detection are done',
        },
        ensembl_annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            id_by => 'ensembl_annotation_build_id',
            is_optional => 1,
        },
        ensembl_annotation_build_id => {
            is => 'Text',
            doc => 'ID of ImportedAnnotation build with the desired ensembl version.',
            default => $ENV{GENOME_DB_ENSEMBL_DEFAULT_IMPORTED_ANNOTATION_BUILD},
        },
    ],
    has_optional_input => [
        dbsnp_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'the dbSNP build by which to determine which variants are novel etc',
        },
        thousand_genomes_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'the 1kg build by which to annotate allele frequencies',
        },
        nhlbi_esp_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'the nhlbi esp build by which to annotate allele frequencies',
        },
        previous_variant_detection_results => {
            is => 'Path',
            doc => 'path to a VCF of previous: skip variant detection and use this',
        },
        nomenclature => {
            is => 'Genome::Nomenclature',
            doc => 'nomenclature used to access clinical data'
        },
        roi_list => {
            is => 'Genome::FeatureList',
            doc => 'only variants in these regions will be included in the final VCF',
        },
        pedigree_file_path => {
            is => 'Genome::File::Tsv',
            doc => 'when supplied overrides the automatic lookup of familial relationships'
        },
        clinical_data_file_path => {
            is => 'Genome::File::Tsv',
            doc => 'when supplied overrides the dumping of clinical data from the nomenclature'
        },
        identify_cases_by => {
            is => 'Text',
            doc => 'the expression which matches "case" samples, typically by their attributes'
        },
        identify_controls_by => {
            is => 'Text',
            doc => 'the expression which matches "control" samples, typically by their attributes'
        },
        #FIXME does this name make sense
        glm_config_file => {
            is => 'Genome::File::Tsv',
            doc => 'the model file for the glm',
        },
    ],
    has_output => [
        multisample_vcf => {
            is => "File",
            is_optional => 1,
            doc => 'Output of variant detection',
        },
        output_directory => {
            is => 'String',
            is_optional => 1,
            doc => 'Output directory for results',
        },
        per_site_report_file => {
            is => "File",
            is_optional => 1,
            doc => "Per site report file generated by joinx vcf-report",
        },
        per_sample_report_file => {
            is => "File",
            is_optional => 1,
            doc => "Per sample report file generated by joinx vcf-report",
        },

    ],
    has_transient_optional => [
        _clinical_data_object => {
            is => "Genome::Model::PhenotypeCorrelation::ClinicalData",
        },
        _sample_names => {
            is => "ARRAY",
            doc => "Sample names found in VCF and population group"
        }
    ],
};

sub help_synopsis_for_create_profile {
    my $self = shift;
    return <<"EOS"

  # quantitative

    genome processing-profile create phenotype-correlation \
      --name 'September 2011 Quantitative Population Phenotype Correlation' \
      --alignment-strategy              'instrument_data aligned to reference_sequence_build using bwa 0.5.9 [-q 5] then merged using picard 1.29 then deduplicated using picard 1.29' \
      --snv-detection-strategy          'samtools r599 filtered by snp-filter v1' \
      --indel-detection-strategy        'samtools r599 filtered by indel-filter v1' \
      --group-samples-for-genotyping-by 'race' \            # some (optional) phenotypic trait, or 'trio' or 'all'
      --phenotype-analysis-strategy     'quantitative' \    # or 'case-control'

    genome propulation-group define 'ASMS-cohort-WUTGI-2011' ASMS1 ASMS2 ASMS3 ASMS4

    genome model define phenotype-correlation \
        --name                      'ASMS-v1' \
        --subject                   'ASMS-cohort-WUTGI-2011' \
        --processing-profile        'September 2011 Quantitative Phenotype Correlation' \


  # case-control

    genome processing-profile create phenotype-correlation \
      --name 'September 2011 Case-Control Population Phenotype Correlation' \
      --alignment-strategy              'instrument_data aligned to reference_sequence_build using bwa 0.5.9 [-q 5] then merged using picard 1.29 then deduplicated using picard 1.29' \
      --snv-detection-strategy          'samtools r599 filtered by snp-filter v1' \
      --indel-detection-strategy        'samtools r599 filtered by indel-filter v1' \
      --group-samples-for-genotyping-by 'trio', \
      --roi_wingspan                    500 \
      --phenotype-analysis-strategy     'case-control'

    genome propulation-group define 'Ceft-Lip-cohort-WUTGI-2011' CL001 CL002 CL003

    genome model define phenotype-correlation \
        --name                  'Cleft-Lip-v1' \
        --subject               'Cleft-Lip-cohort-WUTGI-2011' \
        --processing-profile    'September 2011 Case-Control Phenotype Correlation' \
        --roi_list              'TEST_ROI' \
        --pedigree-file-path    /somedir/somesubdir/thisfamily.ped
        --identify-cases-by     'some_nomenclature.has_cleft_lip = "yes"' \
        --identify-controls-by  'some_nomenclature.has_cleft_lip = "no"' \


    # If you leave off the subject, it would find all patients matching the case/control logic
    # and make a population group called ASMS-v1-cohort automatically???


EOS
}

sub help_detail_for_create_profile {
    return <<EOS
  For a detailed explanation of how to write an alignmen strategy see:
    TBD

  For a detailed explanation of how to write a variant detection strategy, see:
    perldoc Genome::Model::Tools::DetectVariants2::Strategy;

  All builds will have a combined vcf in their variant detection directory.

EOS
}

sub help_manual_for_create_profile {
    return <<EOS
  Manual page content for this pipeline goes here.
EOS
}

#sub _initialize_profile {
sub __profile_errors__ {
    my $class = shift;      # a class method on this model subclass
    my $profile = shift;    # which takes one profile which goes with this sub-type of model and validates it

    my @errors;

    # ensure that each of the variant detection strategies specified will function:
    for my $strategy ('snv','indel','sv','cnv') {
        my $method_name = $strategy . '_detection_strategy';
        if (my $strategy_text = $profile->$method_name) {
            my $strat = Genome::Model::Tools::DetectVariants2::Strategy->get($strategy_text);
            push @errors,
                map {
                    UR::Object::Tag->create(
                        type => 'invalid',
                        properties => [$method_name],
                        desc => $_
                    )
                }
                $strat->__errors__;
        }
    }
    #check joinx
    eval {
        Genome::Model::Tools::Joinx->joinx_path($profile->joinx_version);
    };
    if($@) {
        push @errors, UR::Object::Tag->create( type => 'invalid', properties => [ qw{ joinx_version } ], desc => "Invalid joinx version: $@" );
    }

    return @errors;
}

sub _resolve_resource_requirements_for_build {
    return "-R 'select[mem>8000] rusage[mem=8000]' -M 8000000"
}

our $SHORTCUT_ALIGNMENT_QUERY = 0;
sub _execute_build {
    my ($self,$build) = @_;
    # TODO: remove this and replace with the workflow logic at the bottom when we have one.
    # Version 1 of this pipeline will run in a linear way only if the underlying samples have already
    # had independent alignment and variant detection completed in other models.

    warn "The logic for building this model is only partly functional!  Contact Human Genomics or put in an APIPE-support ticket..";

    #
    # the subject is a population group
    #
    my $population_group = $build->model->subject;
    $build->status_message("subject is " . $population_group->__display_name__);

    #
    # get the reference sequence
    #

    my $reference_build = $build->reference_sequence_build;
    $build->status_message("reference sequence build: " . $reference_build->__display_name__);

    my $annotation_build = $self->ensembl_annotation_build;
    if (!defined $annotation_build) {
        die "No ensembl_annotation_build specified for this model, abort.";
    }
    $build->status_message("ensembl annotation build: " . $annotation_build->__display_name__);

    if (!$annotation_build->is_compatible_with_reference_sequence_build($reference_build)) {
        die "Reference sequence and annotation builds are incompatible!";
    }

    #my $reference_fasta = $reference_build->full_consensus_path('fa');
    #unless(-e $reference_fasta){
    #    die $self->error_message("fasta file for reference build doesn't exist!");
    #}
    #$build->status_message("reference sequence fasta: " . $reference_fasta);

    #
    # get or create the vcf
    #

    $self->output_directory($build->data_directory);


    if ($build->previous_variant_detection_results) {
        #preserve a copy of the multisample VCF in the build directory
        my ($path,$name,$vcf_type) = fileparse($build->previous_variant_detection_results,qr{\.vcf\.gz$},qr{\.vcf$});
        my $copied_filename = $build->data_directory . "/previously_detected_variants" . $vcf_type;
        $self->status_message("Copying build input: " . $build->previous_variant_detection_results. " to $copied_filename");
        my $result = Genome::Sys->copy_file($self->previous_variant_detection_results, $copied_filename);
        unless($result) {
            die $self->error_message("Failed to copy previous variant detection results to the build directory");
        }
        $self->multisample_vcf($copied_filename);
    }
    else {
        my $multisample_vcf = $self->_find_or_generate_multisample_vcf($build,$population_group,$reference_build);
        if($multisample_vcf && -e $multisample_vcf) {
            $self->multisample_vcf($multisample_vcf);

        } else {
            die $self->error_message("$multisample_vcf returned by DV2 but does not exist");
        }
    }
    #notify which VCF
    $self->status_message("Input VCF: " . $self->multisample_vcf);
    $self->status_message("Ensembl annotation build: " . $self->ensembl_annotation_build->name);

    #Do annotation based on available inputs
    my $annotated_vcf = $self->_annotate_multisample_vcf($self->output_directory);
    if($annotated_vcf) {
        $self->status_message("Finished annotating the VCF. The following VCF will be used for further analysis: $annotated_vcf"); #the annotate method updates that
        $self->multisample_vcf($annotated_vcf);
        my $cmd = Genome::Model::Tools::Tabix::Index->create(
            input_file => $annotated_vcf,
            preset => 'vcf',
        );
        if (!$cmd->execute) {
            confess "Failed to create tabix index for file $annotated_vcf";
        }
    }

    # Continue with analysis of the multisample_vcf
    #

    #generate reports

    $self->_generate_callset_reports($self->output_directory);

    #mendelian is broken so let's abort if that's our PP
    if($self->trait_type eq 'mendelian') {
        $self->error_message("Mendelian specific analysis is currently under construction. Model will succeed but not move forwards.");
        return 1;
    }
    if($self->trait_type eq 'none') {
        $self->error_message("No statistical analysis requested.");
        return 1;
    }

    # check that the samples match everywhere they need to
    my $samples = $self->sample_intersection;
    confess "The intersection of samples in the vcf file, population group, and clinical data is empty!" unless @$samples;


    #FIXME Commenting out the delegate class bit until the upstream stuff is sorted out
    my $delegate_class = $self->_generate_delegate_class_name;
    $self->_validate_delegate_class($delegate_class);

    ##For now, simply pass the required items to the delegate class
    my %inputs = $self->_map_properties_to_delegate_inputs($delegate_class);
    my $delegate_obj = $delegate_class->create(%inputs);
    if($delegate_obj) {
        return $delegate_obj->execute;
    }

    return 1;
}

sub _find_or_generate_multisample_vcf {
    my ($self, $build, $population_group, $reference_build) = @_;

    my @builds;
    my @samples;
    # generate or shortcut a multisample VCF

    #
    # get the subject (population group), the individual members and their samples
    #

    my @patients = $population_group->members();
    $build->status_message("found " . scalar(@patients) . " patients");

    @samples = $population_group->samples;
    $build->status_message("found " . scalar(@samples) . " samples");

    my @instdata_assn = $build->inputs(name => 'instrument_data');
    $build->status_message("found " . scalar(@instdata_assn) . " assignments for the current build");

    #my @instdata = Genome::InstrumentData->get(id => [ map { $_->value_id } @instdata_assn ]);
    my @instdata = map { $_->value } @instdata_assn;
    $build->status_message("found " . scalar(@instdata) . " instdata");

    my @per_sample_alignment_results;
    my @bams;

    if ($SHORTCUT_ALIGNMENT_QUERY) {
        # shortcut to speed testing
        @per_sample_alignment_results = Genome::SoftwareResult->get(
            [
            '116553088',
            '116553238',
            '116553281'
            ]
        );

        @builds = Genome::Model::Build->get(
            [
            '116552788',
            '116552996',
            '116553031'
            ]
        );

        @bams = (
            '/gscmnt/gc7001/info/build_merged_alignments/merged-alignment-blade13-4-10.gsc.wustl.edu-rlong-14103-116553088/116553088.bam',
            '/gscmnt/gc7001/info/build_merged_alignments/merged-alignment-blade13-4-10.gsc.wustl.edu-rlong-17210-116553238/116553238.bam',
            '/gscmnt/ams1152/info/build_merged_alignments/merged-alignment-blade13-4-7.gsc.wustl.edu-rlong-12110-116553281/116553281.bam'
        );
    }
    else {
        $self->status_message('Gathering alignments...');
        my $overall_alignment_result = Genome::InstrumentData::Composite->get_or_create(
            inputs => {
                instrument_data => \@instdata,
                reference_sequence_build => $reference_build,
            },
            strategy => $self->alignment_strategy,
            log_directory => $build->log_directory,
        );

        # used by the updated DV2 API
        @per_sample_alignment_results = $overall_alignment_result->_merged_results;
        for my $r (@per_sample_alignment_results) {
            $r->add_user(label => 'uses', user => $build);
        }
        $self->status_message('Found ' . scalar(@per_sample_alignment_results) . ' per-sample alignmnet results.');

        # used directly by the merge tool until we switch to using the above directly
        @bams = $overall_alignment_result->bam_paths;
        $self->status_message('Found ' . scalar(@bams) . ' merged BAMs.');
        for my $bam (@bams){
            unless (-e $bam){
                die $self->error_message("Bam file could not be reached at: ".$bam);
            }
        }

        # this is used by the old, non-DV2 code, but is also used by vcf2maf,
        # which reliese on annotation having been run on the original samples
#        @builds = $self->_get_builds(\@per_sample_alignment_results);

#        my @ar_ids = map { $_->id } @per_sample_alignment_results;
#        my @build_ids = map { $_->id } @builds;
#        print Data::Dumper::Dumper(\@ar_ids, \@build_ids, \@bams);
    }

    #
    # Detect Variants
    #
    # run the DV2 API to do variant detection as we do in somatic, but let it take in N BAMs
    # _internally_ it will (for the first pass):
    #  notice it's running on multiple BAMs
    #  get the single-BAM results
    #  merge them with joinx and make a combined VCF (tolerating the fact that per-bam variants are not VCF)
    #  run bamreadcount to fill-in the blanks
    #

    $self->status_message("Executing detect variants step");

    my %params;
    $params{snv_detection_strategy} = $self->snv_detection_strategy if $self->snv_detection_strategy;
    $params{indel_detection_strategy} = $self->indel_detection_strategy if $self->indel_detection_strategy;
    $params{sv_detection_strategy} = $self->sv_detection_strategy if $self->sv_detection_strategy;
    $params{cnv_detection_strategy} = $self->cnv_detection_strategy if $self->cnv_detection_strategy;

    $params{reference_build_id} = $reference_build->id;

    my $output_dir = $build->data_directory."/variants";
    $params{output_directory} = $output_dir;

    # instead of setting {control_,}aligned_reads_{input,sample}
    # set alignment_results and control_alignment_results

    my @snp_files;
    for my $alignment_result (@per_sample_alignment_results) {
        my $snp_file = $self->get_snp_list_for_something_close_to_this_result($alignment_result);
        if ($snp_file && -e $snp_file) {
            push @snp_files, $snp_file;
        } else {
            $self->error_message("Couldn't locate a variant list from a previous model to speed up QC with for alignment result " . $alignment_result->id . 
                " , (either we no longer make them or the code to find them is broken).");
        }
    }

    unless (scalar (@snp_files) == scalar (@per_sample_alignment_results)) {
        die $self->error_message("Only got " . scalar(@snp_files) . " out of " . scalar(@per_sample_alignment_results) . " expected snplists. Exiting.");
    }




    $params{alignment_results} = \@per_sample_alignment_results;
    $params{control_alignment_results} = [];
    if($build->pedigree_file_path) {
        $params{pedigree_file_path} = $build->pedigree_file_path->path;
        my $ped_is_good_for_polymutt = $self->check_ped_file($build->pedigree_file_path->path);
        if(!$ped_is_good_for_polymutt) {
            $self->error_message("this ped has a missing parent for a person. Unable to proceed with processing, polymutt will crash");
            $self->error_message("Continuing with analysis despite QC failure...");
            $build->add_note(
                header_text => 'QC Failed',
            );
        }
        my $ref_fasta = $reference_build->full_consensus_path("fa");
        my $ped_file = $build->pedigree_file_path->path;
        $self->status_message("About to run Sequencing QC (IBD)");
        $self->status_message("Parameters:");
        $self->status_message("Bams: @bams");
        $self->status_message("Ped file: $ped_file");
        $self->status_message("reference fasta: $ref_fasta");
        $self->status_message("Snp Files: @snp_files");
        $self->status_message("directory: $output_dir/IBD_QC/");
        my $sqc_obj = Genome::Model::Tools::Relationship::SequencingQc->create(
            ped_file=>$ped_file,
            bams=>\@bams,
            reference_fasta=>$reference_build->full_consensus_path("fa"),
            snp_files=>\@snp_files,
            output_dir=>"$output_dir/IBD_QC",
        );
        my $IBD_STATUS = $sqc_obj->execute();
        $sqc_obj->software_result->add_user(label => 'uses', user => $build);
        if($IBD_STATUS eq 'Fail') {
            $self->error_message("Sequencing QC module returned fail code, this ped/model-group has a relationship problem");
            $self->error_message("Continuing with analysis despite QC failure...");
            $build->add_note(
                header_text => 'QC Failed',
            );
        }
        if($IBD_STATUS eq 'Pass') {
            $self->status_message("Sequencing-Qc module reports pass, these individuals seem to be related as the ped has described");
        }

    }
    $params{roi_list} = $build->roi_list;
    $params{roi_wingspan} = $self->roi_wingspan;

    my $command = Genome::Model::Tools::DetectVariants2::Dispatcher->create(%params);
    unless ($command){
        die $self->error_message("Couldn't create detect variants dispatcher from params:\n".Data::Dumper::Dumper \%params);
    }

    my $rv = $command->execute;
    my $err = $@;
    unless ($rv){
        die $self->error_message("Failed to execute detect variants dispatcher(err:$@) with params:\n".Data::Dumper::Dumper \%params);
    }

    $self->status_message("detect variants command completed successfully");

    return $output_dir . '/snvs.vcf.gz';
}

# This method attempts to return a useable snp list (very specific right now, one that ran through samtools->snpfilter) that was previously run on this alignment result.
# If this fails, it falls back to trying to find an alternate alignment result that has a useable snp list.
sub get_snp_list_for_something_close_to_this_result {
    my ($self, $alignment_result) = @_;
    my $return_snp_file;
    my $bam_path = $alignment_result->merged_alignment_bam_path;
    my $reference_build_id = $alignment_result->reference_build->id;

    # Try to find a snplist that was run with this bam
    my @results = Genome::Model::Tools::DetectVariants2::Result::Filter->get(
        filter_name => "Genome::Model::Tools::DetectVariants2::Filter::SnpFilter",
        reference_build_id => $reference_build_id,
        aligned_reads =>$bam_path );
    if(@results) {
        # We are being pretty permissive here. Take ANY result we got.
        my $snp_list = $results[0]->path("snvs.hq.bed");
        $self->status_message("Found a snp list at $snp_list");
        return $snp_list;
    } else {
        $self->status_message("Could not find any DV2::Filter::SnpFilter result object for $bam_path. Widening the search space....");
    }

    # Find a list of other bams that have this same set of instrument data and reference sequence
    my @more_alignment_results = Genome::InstrumentData::AlignmentResult::Merged->get(
        instrument_data_id => [$alignment_result->instrument_data_id],
        reference_build_id => $reference_build_id,
    );

    # For every bam you found, look for an existing useable snp list
    for my $alignment_result (@more_alignment_results) {
        $bam_path = $alignment_result->merged_alignment_bam_path;
        my @results = Genome::Model::Tools::DetectVariants2::Result::Filter->get(
            filter_name => "Genome::Model::Tools::DetectVariants2::Filter::SnpFilter",
            reference_build_id => $reference_build_id,
            aligned_reads =>$bam_path );
        if(@results) {
            # We are being pretty permissive here. Take ANY result we got.
            my $snp_list = $results[0]->path("snvs.hq.bed");
            $self->status_message("Found a snp list at $snp_list");
            return $snp_list;
        }
    }
    return undef;
}


sub _vcf_annotate {
    my ($self, $input_vcf, $annotation_vcf, $info_string) = @_;

    my $output_file = Genome::Sys->create_temp_file_path();
    my $vcf_annotator = Genome::Model::Tools::Joinx::VcfAnnotate->create(
        input_file=> $input_vcf,
        annotation_file=>$annotation_vcf,
        output_file=>$output_file,
        use_bgzip=>1,
        info_fields=>$info_string,
        use_version => $self->joinx_version,
    );
    $vcf_annotator->execute || die "Failed to execute Joinx Vcf annotation using db: $annotation_vcf";
    $self->status_message("Successfully annotated VCF with information from $annotation_vcf");
    return $output_file;
}

sub _dbsnp_info_fields_for_version {
    my $self = shift;
    my $version = shift;
    my @fields = ("GMAF", "dbSNPBuildID=dbSNPBuildID,per-alt", "MUT");
    if ($version >= 137) {
        push(@fields, "PM");
    } else {
        push(@fields, "CLN");
    }
    return join(":", @fields);
}

sub _annotate_multisample_vcf {
    my ($self,$output_dir) = @_;
    my $vcf = $self->multisample_vcf;
    my $annotated_vcf;
    if($self->dbsnp_build) {
        #add dbsnp stuff to the vcf
        #FIXME maybe allow info fields that aren't hardcoded at some point in the future
        my $dbsnp_vcf = $self->dbsnp_build->snvs_vcf;
        $self->status_message("Annotating with dbSNP VCF");
        my $dbsnp_info_fields = $self->_dbsnp_info_fields_for_version($self->dbsnp_build->version);
        $annotated_vcf = $self->_vcf_annotate($vcf, $dbsnp_vcf, $dbsnp_info_fields);
        #set vcf variable so other predefined annotation sources can use it
        $vcf = $annotated_vcf;
    }
    if($self->thousand_genomes_build) {
        #add dbsnp stuff to the vcf
        #FIXME maybe allow info fields that aren't hardcoded at some point in the future
        my $thousand_genomes_vcf = $self->thousand_genomes_build->snvs_vcf;
        $self->status_message("Annotating with 1000 Genomes VCF");
        $annotated_vcf = $self->_vcf_annotate($vcf, $thousand_genomes_vcf, "AF=1KGAF:ASN_AF:AMR_AF:AFR_AF:EUR_AF:SNPSOURCE=1KG,per-alt", $annotated_vcf);
    }
    if($self->nhlbi_esp_build) {
        my $nhlbi_vcf = $self->nhlbi_esp_build->snvs_vcf;
        $self->status_message("Annotating with NHLBI ESP VCF");
        $annotated_vcf = $self->_vcf_annotate($vcf, $nhlbi_vcf, "MAF=NHLBI_ESP_MAF:DP=NHLBI,per-alt", $annotated_vcf);
    }
    if($annotated_vcf) {
        my $new_multisample_vcf = "$output_dir/snvs.merged.annotated.vcf.gz";
        Genome::Sys->copy_file($annotated_vcf, $new_multisample_vcf);
        return $new_multisample_vcf;
    }
    return;
}

sub _generate_callset_reports {
    my ($self, $output_dir) = @_;

    #generate callset reports
    #
    #perl -I ~/src/genome/lib/perl -S gmt joinx vcf-report --info-fields-from-db dbSNPBuildID --input-file snvs.merged.dbSNP135.site_filtered.passing_sites_only.roistrict.vcf.gz --use-bgzip
    my $per_sample_file = "$output_dir/per_sample_callset_report.txt";
    my $per_site_file = "$output_dir/per_site_callset_report.txt";

    #FIXME the input fields should be intelligently set based on available annotation vcfs or taken in from the inputs/params
    #FIXME When vcf report no longer writes reports to Cwd this can go away.
    {
        local $CWD = $output_dir;
        my %params = (
            input_file=> $self->multisample_vcf,
            use_bgzip=>1,
            report_image_type => 'png',
            per_sample_output_file => $per_sample_file,
            per_site_output_file => $per_site_file,
            use_version => $self->joinx_version,
        );
        my @info_fields;
        if($self->dbsnp_build) {
            push @info_fields, 'dbSNPBuildID';
        }
        if($self->thousand_genomes_build) {
            push @info_fields, '1KG';
        }
        if($self->nhlbi_esp_build) {
            push @info_fields, 'NHLBI';
        }
        if(@info_fields) {
            $params{info_fields_from_db} = join(':',@info_fields);
        }
        my $vcf_reporter = Genome::Model::Tools::Joinx::VcfReport->create( %params );
        $vcf_reporter->execute || die "Failed to execute Joinx VCF report";
    }
    #succeeded so set build variables
    $self->per_site_report_file($per_site_file);
    $self->per_sample_report_file($per_sample_file);

    # Compress and index per site report
    my $cmd = Genome::Model::Tools::Tabix::Index->create(
        input_file => $self->_per_site_report_bgzip,
        skip_lines => 1,
        sequence_column => 1,
        start_column => 2,
        end_column => 2,
    );
    if (!$cmd->execute) {
        confess "Failed to create tabix index for file $per_site_file";
    }

    return 1; #not necessary but I feel better knowing it returns true on completion.
}

sub _per_site_report_bgzip {
    my $self = shift;
    my $per_site_file = $self->per_site_report_file;
    my $per_site_file_gz = "$per_site_file.gz";
    return $per_site_file_gz if -s $per_site_file_gz;

    my $bgzip_cmd = "bgzip -c $per_site_file > $per_site_file_gz";
    Genome::Sys->shellcmd(cmd => $bgzip_cmd);
    return $per_site_file_gz;
}

sub _generate_delegate_class_name {
    my ($self) = @_;

    my @trait_words = split('-', $self->trait_type);
    my @cohort_words = split('-', $self->cohort_type);

    my $delegate_class_name = __PACKAGE__ . "::Command::" . join('::', join(q{}, map { ucfirst(lc($_)) } @trait_words), join(q{}, map { ucfirst(lc($_)) } @cohort_words));

    return $delegate_class_name;
}

sub _validate_delegate_class {
    my ($self, $class_name) = @_;
    eval { $class_name->class; };
    if($@) {
        die "Error using class to handle profile with trait " . $self->trait_type . " and cohort type " . $self->cohort_type  . ". $class_name had errors.";
    }
    return 1; #valid class
}

sub _map_properties_to_delegate_inputs {
    my ($self, $class_name) = @_;
    #this maps existing class items to the delegates inputs and reports appropriate warnings
    my @model_inputs = __PACKAGE__->__meta__->properties(is_input => 1);
    my @model_input_names = map { $_->property_name } @model_inputs;
    my %build_input_hash = map { $_ => $self->$_ } @model_input_names;

    #add inputs
    #right now this is the multisample snv file
    #$build_input_hash{multisample_vcf} = $self->_multisample_vcf;

    #get delegate class inputs
    my @delegate_inputs = $class_name->__meta__->properties(is_input => 1);
    my @delegate_input_names = map { $_->property_name } @delegate_inputs;
    my %delegate_input_hash;
    #assuming all are REQUIRED for now. FIXME later
    for my $expected_input (@delegate_input_names) {
        if(defined($build_input_hash{$expected_input})) {
            #all is well
            $delegate_input_hash{$expected_input} = delete $build_input_hash{$expected_input};
        }
        else {
            #input is either undefined OR not listed
            if($self->can($expected_input)) {
                # Assume this is either needing to be run or already set as an output variable
                # only assign it if defined
                my $value = $self->$expected_input;
                $delegate_input_hash{$expected_input} = $value if defined $value;
            }
            else {
                die "Expected input $expected_input of delegate class $class_name not available. Check that this is specified as a model input to your model."
            }
        }
    }

    #report excess inputs
    for my $available_input (keys %build_input_hash) {
        $self->error_message("$available_input unused by $class_name and will be ignored.");
    }
    return %delegate_input_hash;
}

sub clinical_data_object {
    my $self = shift;

    return $self->_clinical_data_object if $self->_clinical_data_object;

    if ($self->clinical_data_file_path) {
        $self->status_message("Loading clinical data from file " . $self->clinical_data_file_path->path);

        $self->_clinical_data_object(Genome::Model::PhenotypeCorrelation::ClinicalData->from_file(
            $self->clinical_data_file_path->path
            ));
    } else {
        $self->status_message("Loading clinical data from database using nomenclature " . $self->nomenclature->name);
        $self->_clinical_data_object(Genome::Model::PhenotypeCorrelation::ClinicalData->from_database(
            $self->nomenclature,
            $self->subject->samples
            ));
    }

    return $self->_clinical_data_object;
}

sub clinical_data_file {
    my ($self) = @_;
    my $cd = $self->clinical_data_object();
    $self->status_message("preparing clinical data files\n");

    my $clinical_data = $self->output_directory . "/Clinical_Data.txt";
    my $clinical_data_md5 = $self->output_directory . "/Clinical_Data.txt.md5";

    $self->status_message("Writing clinical data to $clinical_data");
    my $digest = $cd->to_file($clinical_data);
    $self->status_message("m5sum of input clinical data: " . $digest . "\n");
    my $fh = Genome::Sys->open_file_for_writing($clinical_data_md5);
    $fh->write("$digest\n");
    $fh->close;
    return $clinical_data;
}

sub sample_intersection {
    my $self = shift;
    if (!$self->_sample_names) {
        my %pop_group_samples = map {$_->name => 1} $self->subject->samples;
        my @vcf_samples = $self->_samples_from_vcf;
        my @clinical_samples = $self->clinical_data_object->sample_names;
        $self->status_message("Population group contains " . scalar keys(%pop_group_samples) . " samples");
        $self->status_message("Multisample vcf contains " . scalar @vcf_samples . " samples");
        $self->status_message("Clinical data file contains " . scalar @clinical_samples . " samples");

        # compute the intersection of samples in population group, vcf, and clinical data file
        my $samples_hash = reduce {
            our ($a, $b);
            return { map {$_ => 1} grep {exists $a->{$_}} @$b }
            } \%pop_group_samples, \@vcf_samples, \@clinical_samples;

        my @samples = sort keys %$samples_hash;
        $self->status_message("The sample intersection contains " . scalar @samples . " samples");

        $self->_sample_names(\@samples);
    }
    return $self->_sample_names;
}

sub sample_list_file {
    my ($self) = @_;

    my $sample_file = $self->output_directory . "/Sample_List.txt";
    $self->status_message("Attempting to generate sample list file: $sample_file");
    my $sample_fh = Genome::Sys->open_file_for_writing($sample_file);

    my $sample_names = $self->sample_intersection;
    print $sample_fh join("\n",@$sample_names),"\n";
    close($sample_fh);
    return $sample_file;
}

sub glm_model_file {
    my ($self) = @_;
    $self->status_message("preparing glm model file\n");
    my $glm_model = $self->output_directory . "/glm-model.txt";
    if(defined($self->glm_config_file)) {
        #copy it over
        $self->status_message("Copying build input: " . $self->glm_config_file->path . " to $glm_model");
        Genome::Sys->copy_file($self->glm_config_file->path,$glm_model);    #this croaks if it fails
        return $glm_model;
    }
    else {
        die "Autogeneration of glm file currently unimplemented";
    }
}

#FIXME It would be better to either add this to the PP
sub glm_max_cols_per_file {
    my ($self) = @_;
    return 1000;
}

sub _samples_from_vcf {
    my ($self) = @_;
    my $reader = Genome::File::Vcf::Reader->new($self->multisample_vcf);
    return $reader->header->sample_names;
}

sub _get_builds {
    my $self = shift;
    my $results = shift;
    my @results = @{$results};
    my @builds;
    for my $result (@results) {
        my @users_who_are_builds = grep { $_->user_class_name =~ m/Genome\:\:Model\:\:Build\:\:ReferenceAlignment/ } $result->users;
        push @builds, Genome::Model::Build->get($users_who_are_builds[0]->user_id);
    }
    return @builds;
}


sub _validate_build {
    # this is where we sanity check things like inputs making sense before actually building
    my $self = shift;
    my $dir = $self->data_directory;

    my @errors;
    unless (1) {
        my $e = $self->error_message("Something is wrong!");
        push @errors, $e;
    }

    if (@errors) {
        return;
    }
    else {
        return 1;
    }
}


sub check_ped_file {
    my ($self, $ped_file) = @_;
    my $fh = IO::File->new($ped_file);
    my $complete_trio=0;
    while(my $line = $fh->getline) {
        my ($family_id, $ind_id, $dad, $mom, $sex, $glf, $affect)  = split "\t", $line;
        if($dad && !$mom) {
            return 0;
        }
        if(!$dad && $mom) {
            return 0;
        }
        if($dad && $mom) {
            $complete_trio=1;
        }
    }
    return $complete_trio;
}



1;

__END__

# TODO: replace the above _execute_build with an actual workflow

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

    #### This is old code from the somatic variation pipeline, replace with phenotype correlation params/inputs! #####

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

    push @inputs, build_id => $build->id;

    return @inputs;
}


1;
