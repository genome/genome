package Genome::Model::GenotypeMicroarray;

use strict;
use warnings;

use Genome;
use File::Basename qw(dirname basename);

class Genome::Model::GenotypeMicroarray{
    is => 'Genome::ModelDeprecated',
    has_input => {
        dbsnp_build_id => {
            is => 'Text',
            to => 'value_id',
            where => [ name => 'dbsnp_build', value_class_name => 'Genome::Model::Build::ImportedVariationList' ],
            is_many => 0,
            is_mutable => 1,
            is_optional => 1,
            doc => 'dbsnp build that this model is built against'
        },
    },
    has => {
        dbsnp_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            id_by => 'dbsnp_build_id',
        },
        dbsnp_version => { 
            is => 'Text',
            via => 'dbsnp_build',
            to => 'version',
        },
        reference_sequence_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            via => 'dbsnp_build', 
            to => 'reference'
        },
        reference_sequence_build_id => {
            is => 'Text',
            via => 'reference_sequence_build',
            to => 'id',
        },
        refseq_name => { 
            is => 'Text',
            via => 'reference_sequence_build',
            to => 'name',
        },
        refseq_version => { 
            is => 'Text',
            via => 'reference_sequence_build',
            to => 'version',
        },
    },
    has_calculated => {
        genotype_vcf => { is => 'Text', doc => 'Returns the VCF file name if exists, or NA if not.', },
    },
};

sub sequencing_platform { return 'genotype file'; }

sub is_internal { 
    my $self = shift;
    my ($instrument_data) = $self->instrument_data;
    my $source = $instrument_data->import_source_name;
    if (defined $source and $source =~ /wugc/i) {
        return 1;
    }
    return 0;
}

sub _additional_parts_for_default_name {
    my ($self, %params) = @_;
    my ($instrument_data) = $self->instrument_data;
    if ( not $instrument_data ) {
        $instrument_data = $params{instrument_data};
        if ( not $instrument_data ) {
            die 'No instrument data found for model';
        }
    }
    return ( $instrument_data->import_source_name, $instrument_data->sequencing_platform, $self->refseq_name );
}

my @format_types = (
    {
        id => 'GT',
        name => 'genotype',
        header => ',Number=1,Type=String,Description="Genotype">',
    },
    {
        id => 'ALLELES',
        name => 'alleles',
        header => ',Number=1,Type=String,Description="Alleles called from the microarray chip">',
    },
    {
        id => 'BAF',
        name => 'b_allele_freq',
        header => ',Number=1,Type=Float,Description="B allele frequency">',
    },
    {
        id => 'CNV_CONF',
        name => 'cnv_confidence',
        header => ',Number=1,Type=Float,Description="CNV Confidence">',
    },
    {
        id => 'CNV_VAL',
        name => 'cnv_value',
        header => ',Number=1,Type=Float,Description="CNV Value">',
    },
    {
        id => 'LOG_R',
        name => 'log_r_ratio',
        header => ',Number=1,Type=Float,Description="Log R Ratio">',
    },
    {
        id => 'GC_SCORE',
        name => 'gc_score',
        header => ',Number=1,Type=Float,Description="GC Score">',
    },
);

sub format_types {
    return @format_types;
}

sub format_name_for_id {
    my ($self, $format_type_id) = @_;
    my $format = List::Util::first { $_->{id} eq $format_type_id }  @format_types;
    return if not $format;
    return $format->{name};
}

sub genotype_filters {
    # FIXME put in processing profile [or the like]
    my $self = shift;

    my @filters = (qw/ gc_score:min=0.7 /); 
    push @filters, 'invalid_iscan_ids' if $self->reference_sequence_build->version eq '36';

    return \@filters;
}

sub genotype_vcf {
    my $self = shift;

    my $build = $self->last_succeeded_build;
    return 'NO_SUCCEESSFUL_BUILD' if not $build;

    my $vcf = $build->original_genotype_vcf_file_path;
    return 'NO_VCF_FOR_BUILD' if not -s $vcf;

    return $vcf;
}

sub dependent_cron_ref_align {
    my $self = shift;

    my @subjects = ($self->subject);
    push @subjects, Genome::Sample->get(default_genotype_data_id => [map { $_->id } $self->instrument_data]);

    my @ref_align_models = Genome::Model::ReferenceAlignment->get(
        subject_id => [map { $_->id } @subjects],
    );

    # limit to models with a compatible reference sequence build
    my $gm_rsb = $self->reference_sequence_build;
    my @compatible_ref_align_models = grep {
        my $ra_rsb = $_->reference_sequence_build;
        $ra_rsb->is_compatible_with($gm_rsb);
    } @ref_align_models;

    # limit to models that either don't have a genotype_microarray_model yet or have the same genotype_microarray_model
    my @dependent_models = grep {
        my $gmm = $_->genotype_microarray_model;
        (!$gmm or ($gmm && $gmm->id eq $self->id));
    } @compatible_ref_align_models;

    # limit to models for which new results might still be useful
    my @current_models = grep {
        $_->config_profile_item and $_->config_profile_item->is_current
    } @dependent_models;

    return @current_models;
}

sub request_builds_for_dependent_cron_ref_align {
    my $self = shift;
    my $sample = $self->subject;
    return 1 unless $sample->class eq 'Genome::Sample';

    for my $ref_align ($self->dependent_cron_ref_align) {
        my @lane_qc = $ref_align->get_lane_qc_models;
        for (@lane_qc) { $_->build_requested(1) };
        $ref_align->build_requested(1);
    }
    return 1;
}

sub _resolve_resource_requirements_for_build {
    return "-R 'select[mem>4000] rusage[mem=4000]' -M 4000000"
}

#< Work Flow >#
sub map_workflow_inputs {
    my ($self, $build) = @_;
    my @instrument_data = $build->instrument_data;
    return (
        build => $build,
    );
}

sub _resolve_workflow_for_build {
    my ($self, $build, $lsf_queue, $lsf_project) = @_;

    $lsf_queue = Genome::Config::get('lsf_queue_build_worker_alt') unless defined($lsf_queue);
    $lsf_project = 'build' . $build->id unless defined($lsf_project);

    my $workflow = Genome::WorkflowBuilder::DAG->create(
        name => $build->workflow_name,
    );

    my $previous_op;
    my $add_operation = sub{
        my ($name) = @_;
        my $command_class_name = 'Genome::Model::GenotypeMicroarray::Build::'.join('', map { ucfirst } split(' ', $name));

        my $operation = Genome::WorkflowBuilder::Command->create(
            name => $name,
            command => $command_class_name,
        );
        $operation->lsf_queue($lsf_queue);
        $operation->lsf_project($lsf_project);

        $workflow->add_operation($operation);

        if ($previous_op) {
            $workflow->create_link(
                source => $previous_op,
                source_property => 'build',
                destination => $operation,
                destination_property => 'build',
            );
        }

        return $operation;
    };

    my $create_og_files_op = $add_operation->('create original genotype files');
    $previous_op = $create_og_files_op;

    $workflow->connect_input(
        input_property => 'build',
        destination => $create_og_files_op,
        destination_property => 'build',
    );

    my $create_filtered_genotype_file_op = $add_operation->('create filtered genotype tsv file');
    $previous_op = $create_filtered_genotype_file_op;

    my $create_copy_number_file = $add_operation->('create copy number tsv file');
    $previous_op = $create_copy_number_file;

    my $create_gold_snp_file_op = $add_operation->('create gold snp file');
    $previous_op = $create_gold_snp_file_op;

    my $create_gold_snp_bed_file_op = $add_operation->('create gold snp bed file');
    $previous_op = $create_gold_snp_bed_file_op;

    $workflow->connect_output(
        source => $previous_op,
        source_property => 'build',
        output_property => 'build',
    );

    $workflow->recursively_set_log_dir($build->log_directory);

    return $workflow;
}
#<>#

1;

