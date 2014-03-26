package Genome::Model::ClinSeq::TestData;

use strict;
use warnings;
use Genome;

use Genome::Test::Factory::Build;
use Genome::Test::Factory::Library;
use Genome::Test::Factory::Model::ClinSeq;
use Genome::Test::Factory::Model::DifferentialExpression;
use Genome::Test::Factory::Model::ImportedAnnotation;
use Genome::Test::Factory::Model::ImportedReferenceSequence;
use Genome::Test::Factory::Model::ImportedVariationList;
use Genome::Test::Factory::Model::ReferenceAlignment;
use Genome::Test::Factory::Model::RnaSeq;
use Genome::Test::Factory::Model::SomaticVariation;
use Genome::Test::Factory::ProcessingProfile::ClinSeq;
use Genome::Test::Factory::ProcessingProfile::DifferentialExpression;
use Genome::Test::Factory::ProcessingProfile::ReferenceAlignment;
use Genome::Test::Factory::ProcessingProfile::RnaSeq;
use Genome::Test::Factory::ProcessingProfile::SomaticVariation;
use Genome::Utility::Test;

sub load {
    my $self = shift;
    my %params = @_;
    my %ids;
    #default base_dir
    my $base_dir = $ENV{GENOME_TEST_INPUTS}."Genome-Model-ClinSeq-TestData/2014-03-26";

    if ($params{base_dir}) {
        $base_dir = delete $params{base_dir};
    }

    $ENV{GENOME_DB} = "$base_dir/reference_annotations/";

    my $individual = Genome::Individual->create(common_name => "FAKE1",
        name => "test-clin-seq",
        gender => "unspecified",
        upn => "-353",
    );
    $ids{TEST_INDIVIDUAL_ID} = $individual->id;

#Obtain a normal DNA sample
    my $normal_dna_sample = Genome::Sample->create(source => $individual,
        name => "clinseq-normal-dna",
        common_name => "normal",
        extraction_type => "genomic dna",
        cell_type => "primary",
        tissue_desc => "blood",
    );
    my $normal_dna_sample_id = $normal_dna_sample->id;
    $ids{NORMAL_DNA_SAMPLE} = $normal_dna_sample_id;
    my $normal_inst_data = create_instrument_data_from_sample($normal_dna_sample);

#Obtain a normal RNA sample
    my $normal_rna_sample = Genome::Sample->create(source => $individual,
        name => "clinseq_normal_rna",
        common_name => "normal",
        extraction_type => "rna",
        cell_type => "primary",
        tissue_desc => "brain",
    );
    my $normal_rna_sample_id = $normal_rna_sample->id;
    $ids{NORMAL_RNA_SAMPLE} = $normal_rna_sample_id;
    my $normal_rna_inst_data = create_instrument_data_from_sample($normal_rna_sample);


#Obtain a tumor DNA sample
    my $tumor_dna_sample = Genome::Sample->create(source => $individual,
        name => "clinseq-tumor-dna",
        common_name => "met",
        extraction_type => "genomic dna",
        cell_type => "primary",
        tissue_desc => "brain",
    );
    my $tumor_dna_sample_id = $tumor_dna_sample->id;
    $ids{TUMOR_DNA_SAMPLE} = $tumor_dna_sample_id;
    my $tumor_inst_data = create_instrument_data_from_sample($tumor_dna_sample);

#Obtain a tumor RNA sample
    my $tumor_rna_sample = Genome::Sample->create(source => $individual,
        name => "clinseq_tumor_rna",
        common_name => "met",
        extraction_type => "rna",
        cell_type => "primary",
        tissue_desc => "brain",
    );
    my $tumor_rna_sample_id = $tumor_rna_sample->id;
    $ids{TUMOR_RNA_SAMPLE} = $tumor_rna_sample_id;
    my $rna_inst_data = create_instrument_data_from_sample($tumor_rna_sample);

    my $dbsnp_model = Genome::Test::Factory::Model::ImportedVariationList->setup_object;
    my $dbsnp_build = Genome::Test::Factory::Build->setup_object(model_id => $dbsnp_model->id);
    $ids{DBSNP_BUILD} = $dbsnp_build->id;

    my $annotation_model = Genome::Test::Factory::Model::ImportedAnnotation->setup_object;
    my $annotation_build = Genome::Test::Factory::Build->setup_object(model_id => $annotation_model->id,
                                                                         data_directory => $base_dir."/annot_dir",
                                                                         status => "Succeeded",
                                                                        );
    $ids{ANNOTATION_BUILD} = $annotation_build->id;

    my $ref_seq_model = Genome::Test::Factory::Model::ImportedReferenceSequence->setup_object;
    my $ref_seq_build = Genome::Test::Factory::Build->setup_object(model_id => $ref_seq_model->id);
    $ids{REFSEQ_BUILD} = $ref_seq_build->id;

#setup refalign build
    my $ref_align_pp = Genome::Test::Factory::ProcessingProfile::ReferenceAlignment->setup_object();
    $ids{REFALIGN_PP} = $ref_align_pp->id;
    my $normal_model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(
        reference_sequence_build => $ref_seq_build,
        dbsnp_build => $dbsnp_build,
        annotation_reference_build => $annotation_build,
        processing_profile_id => $ref_align_pp->id,
        subject_id => $normal_dna_sample->id,
    );
    $normal_model->add_instrument_data($normal_inst_data);
    $ids{NORMAL_REFALIGN_MODEL} = $normal_model->id;
    my $normal_build = Genome::Test::Factory::Build->setup_object(model_id => $normal_model->id, status => "Succeeded");
    my $tumor_model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(
        reference_sequence_build => $ref_seq_build,
        dbsnp_build => $dbsnp_build,
        annotation_reference_build => $annotation_build,
        processing_profile_id => $ref_align_pp->id,
        subject_id => $tumor_dna_sample->id,
    );
    $tumor_model->add_instrument_data($tumor_inst_data);
    $ids{TUMOR_REFALIGN_MODEL} = $tumor_model->id;
    my $tumor_build = Genome::Test::Factory::Build->setup_object(model_id => $tumor_model->id, status => 'Succeeded');
    my %clinseq_model_params;
    unless ($params{exclude_wgs_model}) {
        my $wgs_pp = Genome::Test::Factory::ProcessingProfile::SomaticVariation->setup_object();
        $ids{WGS_PP} = $wgs_pp->id;
        my $wgs_model = Genome::Test::Factory::Model::SomaticVariation->setup_object(
            subject_id => $tumor_model->subject->id,
            normal_model => $normal_model,
            tumor_model => $tumor_model,
            processing_profile_id => $wgs_pp->id,
            annotation_build => $annotation_build,
            previously_discovered_variations_build => $dbsnp_build,
        );
        $ids{WGS_MODEL} = $wgs_model->id;
        my $wgs_build = Genome::Test::Factory::Build->setup_object(model_id => $wgs_model->id, status => 'Succeeded',
            data_directory => $base_dir."/som_var_dir",
        );
        $ids{WGS_BUILD} = $wgs_build->id;
        $clinseq_model_params{wgs_model} = $wgs_model;
    }
    my $exome_pp = Genome::Test::Factory::ProcessingProfile::SomaticVariation->setup_object(
        snv_detection_strategy => "strelka 0.4.6.2 [isSkipDepthFilters = 0]");
    $ids{EXOME_PP} = $exome_pp->id;

    unless ($params{exclude_exome_model}) {
            my $exome_model = Genome::Test::Factory::Model::SomaticVariation->setup_object(
            subject_id => $tumor_model->subject->id,
            normal_model => $normal_model,
            tumor_model => $tumor_model,
            processing_profile_id => $exome_pp->id,
            annotation_build => $annotation_build,
            previously_discovered_variations_build => $dbsnp_build,
        );
        $ids{EXOME_MODEL} = $exome_model->id;
        my $exome_build = Genome::Test::Factory::Build->setup_object(model_id => $exome_model->id, status => 'Succeeded',
            data_directory => $base_dir."/exome_som_var_dir",
        );
        $ids{EXOME_BUILD} = $exome_build->id;
        $clinseq_model_params{exome_model} = $exome_model;
    }

    my $rna_seq_pp = Genome::Test::Factory::ProcessingProfile::RnaSeq->setup_object();
    $ids{RNASEQ_PP} = $rna_seq_pp->id;
    
    my $cancer_annotation_db_id = Genome::Model::ClinSeq->__meta__->property("cancer_annotation_db")->default_value;

    unless ($params{exclude_normal_rnaseq_model}) {
        my $rna_seq_model = Genome::Test::Factory::Model::RnaSeq->setup_object(
            subject_id => $normal_rna_sample->id,
            processing_profile_id => $rna_seq_pp->id,
            annotation_build => $annotation_build,
            reference_sequence_build => $ref_seq_build,
            cancer_annotation_db => $cancer_annotation_db_id,
        );
        $rna_seq_model->add_instrument_data($normal_rna_inst_data);
        my $rna_seq_build = Genome::Test::Factory::Build->setup_object(model_id => $rna_seq_model->id, 
                                                                        status => "Succeeded",
                                                                        data_directory => "$base_dir/rnaseq_dir");
        $ids{RNASEQ_MODEL} = $rna_seq_model->id;
        $clinseq_model_params{normal_rnaseq_model} = $rna_seq_model;
    }
    unless ($params{exclude_tumor_rna_seq_model}) {
        my $rna_seq_model = Genome::Test::Factory::Model::RnaSeq->setup_object(
            subject_id => $tumor_rna_sample->id,
            processing_profile_id => $rna_seq_pp->id,
            annotation_build => $annotation_build,
            reference_sequence_build => $ref_seq_build,
            cancer_annotation_db => $cancer_annotation_db_id,
        );
        $rna_seq_model->add_instrument_data($rna_inst_data);
        $ids{TUMOR_RNASEQ_MODEL} = $rna_seq_model->id;
        my $rna_seq_build = Genome::Test::Factory::Build->setup_object(model_id => $ids{TUMOR_RNASEQ_MODEL}, 
                                                                        status => 'Succeeded',
                                                                        data_directory => "$base_dir/tumor_rnaseq_dir");
        $clinseq_model_params{tumor_rnaseq_model} = $rna_seq_model;
    }
    my $diff_ex_pp = Genome::Test::Factory::ProcessingProfile::DifferentialExpression->setup_object;
    $ids{DIFFEXP_PP} = $diff_ex_pp->id;
    unless ($params{exclude_de_model}) {
        my $diff_exp_model = Genome::Test::Factory::Model::DifferentialExpression->setup_differential_expression_model(
            processing_profile_id => $ids{DIFFEXP_PP},
            normal_model_id => $ids{RNASEQ_MODEL},
            tumor_model_id => $ids{TUMOR_RNASEQ_MODEL},
        );
        $ids{DIFFEXP_MODEL} = $diff_exp_model->id;
        my $diff_exp_build = Genome::Test::Factory::Build->setup_object(
            model_id => $ids{DIFFEXP_MODEL},
            status => "Succeeded",
            data_directory => "$base_dir/de_dir",
        );
        $ids{DIFFEXP_BUILD} = $diff_exp_build->id;
        $clinseq_model_params{de_model} = $diff_exp_model;
    }
    my $clin_seq_pp = Genome::Test::Factory::ProcessingProfile::ClinSeq->setup_object;
    $ids{CLINSEQ_PP} = $clin_seq_pp->id;
    
    my $clin_seq_model = Genome::Test::Factory::Model::ClinSeq->setup_object(
        %clinseq_model_params,
        processing_profile_id => $clin_seq_pp->id,
    );
    $ids{CLINSEQ_MODEL} = $clin_seq_model->id;
    my $clin_seq_build = Genome::Test::Factory::Build->setup_object(model_id => $clin_seq_model->id, status => "Succeeded", data_directory => "$base_dir/clinseq_dir");
    $ids{CLINSEQ_BUILD} = $clin_seq_build->id;
    $ids{CANCER_ANNOT_PATH} = $clin_seq_build->cancer_annotation_db->data_directory;
    $ids{MISC_ANNOT_PATH} = $clin_seq_build->misc_annotation_db->data_directory;
    $ids{COSMIC_ANNOT_PATH} = $clin_seq_build->cosmic_annotation_db->data_directory;
    return \%ids;
}

sub create_instrument_data_from_sample {
    my $sample = shift;
    my $lib = Genome::Test::Factory::Library->setup_object(
        sample_id => $sample->id,
    );
    my $inst_data = Genome::InstrumentData::Solexa->create(library => $lib);
    return $inst_data;
}

1;

