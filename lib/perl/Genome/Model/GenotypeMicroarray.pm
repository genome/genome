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
    has_param => {
        input_format => {
            doc => 'file format, defaults to "wugc", which is currently the only format supported',
            valid_values => ['wugc'],
            default_value => 'wugc',
        },
        instrument_type => {
            doc => 'the type of microarray instrument',
            valid_values => [qw/ affymetrix illumina infinium plink unknown /],
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

our $format_types = {
    GT => {
        id => 'GT',
        name => 'genotype',
        header => ',Number=1,Type=String,Description="Genotype">',
    },
    ALLELES => {
        id => 'ALLELES',
        name => 'alleles',
        header => ',Number=1,Type=String,Description="Alleles called from the microarray chip">',
    },
    CNV_CONF => {
        id => 'CNV_CONF',
        name => 'cnv_confidence',
        header => ',Number=1,Type=Float,Description="CNV Confidence">',
    },
    CNV_VAL => {
        id => 'CNV_VAL',
        name => 'cnv_value',
        header => ',Number=1,Type=Float,Description="CNV Value">',
    },
    LOG_R => {
        id => 'LOG_R',
        name => 'log_r_ratio',
        header => ',Number=1,Type=Float,Description="Log R Ratio">',
    },
    GC_SCORE => {
        id => 'GC_SCORE',
        name => 'gc_score',
        header => ',Number=1,Type=Float,Description="GC Score">',
    },
};
sub format_types {
    return values %$format_types;
}

sub format_name_for_id {
    return $format_types->{$_[1]}->{name};
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
        reference_sequence_build => $self->reference_sequence_build,
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
        (not $gmm || ($gmm && $gmm->id == $self->id));
    } @compatible_ref_align_models;

    return @dependent_models;
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

sub _execute_build {
    my ($self, $build) = @_;
    $self->debug_message('Execute genotype microarray build '.$build->__display_name__);

    my $instrument_data = $build->instrument_data;
    if ( not $instrument_data ) {
        $self->error_message('No instrument data for genotype microarray build '.$build->__display_name__);
        return;
    }
    $self->debug_message('Instrument data: '.$instrument_data->id.' '.$instrument_data->sequencing_platform);

    my $reference_sequence_build = $build->model->reference_sequence_build;
    if ( not $reference_sequence_build ) {
        $self->error_message('No reference sequence build for '.$build->__display_name__);
        return;
    }
    $self->debug_message('Reference sequence build: '.$reference_sequence_build->__display_name__);

    my $dbsnp_build = $build->dbsnp_build;
    if ( not $dbsnp_build ) {
        $dbsnp_build = Genome::Model::ImportedVariationList->dbsnp_build_for_reference($reference_sequence_build);
        if ( not $dbsnp_build ) {
            $self->error_message('No dbsnp build for '.$build->__display_name__);
            return;
        }
        $build->dbsnp_build($dbsnp_build);
        $build->model->dbsnp_build($dbsnp_build);
    }
    $self->debug_message('DB SNP build: '.$dbsnp_build->__display_name__);

    ###
    # Original genotype files: VCF and TSV
    my $create_og_files = Genome::Model::GenotypeMicroarray::Build::CreateOriginalGenotypeFiles->create(
        build => $build,
    );
    if ( not $create_og_files ) {
        $self->error_message('Failed to create command to create extract command original genotype VCF file!');
        return;
    }
    $create_og_files->dump_status_messages(1);
    if ( not $create_og_files->execute ) {
        $self->error_message('Failed to execute command to create extract command original genotype VCF file!');
        return;
    }

    ###
    # Filters for extracting from the above original file
    my @filters = (qw/ gc_score:min=0.7 /); 
    push @filters, 'invalid_iscan_ids' if $reference_sequence_build->version eq '36';

    my $create_filtered_genotype_file = Genome::Model::GenotypeMicroarray::Build::CreateFilteredGenotypeTsvFile->create(
        build => $build,
    );
    if ( not $create_filtered_genotype_file ) {
        $self->error_message('Failed to create command to create genotype file!');
        return;
    }
    $create_filtered_genotype_file->dump_status_messages(1);
    if ( not $create_filtered_genotype_file->execute ) {
        $self->error_message('Failed to execute command to create genotype file!');
        return;
    }

    # Copy number file. No headers, tab sep with chrom, pos and log r ratio
    my $create_copy_number_file = Genome::Model::GenotypeMicroarray::Build::CreateCopyNumberTsvFile->create(
        build => $build,
    );
    if ( not $create_copy_number_file ) {
        $self->error_message('Failed to create command to create copy number file!');
        return;
    }
    $create_copy_number_file->dump_status_messages(1);
    if ( not $create_copy_number_file->execute ) {
        $self->error_message('Failed to execute command to create copy number file!');
        return;
    }

    my $gold_snp_cmd = Genome::Model::GenotypeMicroarray::Build::CreateGoldSnpFile->create(
        build => $build,
    );
    if ( not $gold_snp_cmd ) {
        $self->error_message("Cannot create gold snp tool.");
        return;
    }
    $gold_snp_cmd->dump_status_messages(1);
    if ( not $gold_snp_cmd->execute ) {
        $self->error_message("Cannot execute gold snp tool");
        return;
    }

    $self->debug_message('Create gold snp bed file...');
    my $snvs_bed = $build->snvs_bed;
    $self->debug_message('Gold snp bed file: '.$snvs_bed);
    my $gold_snp_bed = Genome::Model::GenotypeMicroarray::Command::CreateGoldSnpBed->create(
        input_file => $build->formatted_genotype_file_path,
        output_file => $snvs_bed,
        reference => $reference_sequence_build,
    );
    if ( not $gold_snp_bed ) {
        $self->error_message('Failed to create gold snp bed tool!');
        return;
    }
    $gold_snp_bed->dump_status_messages(1);
    unless ($gold_snp_bed->execute) {
        $self->error_message("Could not generate gold snp bed file at $snvs_bed from snp array file");
        return;
    }
    if ( not -s $snvs_bed ) {
        $self->error_message("Executed 'create gold snp bed', but snvs bed file ($snvs_bed) does not exist");
        return;
    }
    $self->debug_message("Create gold snp bed file...OK");

    $self->debug_message('Execute genotype microarray build...OK');
    return 1;
}

1;

