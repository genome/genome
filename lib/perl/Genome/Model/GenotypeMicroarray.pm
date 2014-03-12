package Genome::Model::GenotypeMicroarray;

use strict;
use warnings;

use Genome;
use File::Basename qw(dirname basename);

class Genome::Model::GenotypeMicroarray{
    is => 'Genome::ModelDeprecated',
    has_param => [
        input_format => {
            doc => 'file format, defaults to "wugc", which is currently the only format supported',
            valid_values => ['wugc'],
            default_value => 'wugc',
        },
        instrument_type => {
            doc => 'the type of microarray instrument',
            valid_values => [qw/ affymetrix illumina infinium plink unknown /],
        },
    ],
    has => [
        reference_sequence_build_id => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'reference_sequence_build', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence' ],
            is_many => 0,
            is_mutable => 1, # TODO: make this non-optional once backfilling is complete and reference placeholder is deleted
            is_optional => 1,
            doc => 'reference sequence to align against'
        },
        reference_sequence_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            id_by => 'reference_sequence_build_id',
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
        dbsnp_build_id => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'dbsnp_build', value_class_name => 'Genome::Model::Build::ImportedVariationList' ],
            is_many => 0,
            is_mutable => 1,
            is_optional => 1,
            doc => 'dbsnp build that this model is built against'
        },
        dbsnp_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            id_by => 'dbsnp_build_id',
        },
        dbsnp_version => { 
            is => 'Text',
            via => 'dbsnp_build',
            to => 'version',
        },
    ],
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

sub dependent_cron_ref_align {
    my $self = shift;

    my @subjects = ($self->subject);
    push @subjects, Genome::Sample->get(default_genotype_data_id => [map { $_->id } $self->instrument_data]);

    my @ref_align_models = Genome::Model::ReferenceAlignment->get(
        subject_id => [map { $_->id } @subjects],
        reference_sequence_build => $self->reference_sequence_build,
        auto_assign_inst_data => 1, # our current way of saying auto-build, later to be a project relationship
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

    # Original genotype VCF file
    $self->debug_message('Create original genotype VCF file...');
    my $original_genotype_vcf_file = $build->original_genotype_vcf_file_path;
    $self->debug_message('Original genotype file: '.$original_genotype_vcf_file);
    my $extract = Genome::Model::GenotypeMicroarray::Command::Extract->create(
        instrument_data => $instrument_data,
        variation_list_build => $dbsnp_build,
        output => $original_genotype_vcf_file,
    );
    if ( not $extract ) {
        $self->error_message('Failed to create command to create extract command original genotype VCF file!');
        return;
    }
    $extract->dump_status_messages(1);
    if ( not $extract->execute ) {
        $self->error_message('Failed to execute command to create extract command original genotype VCF file!');
        return;
    }
    # Check that genotypes were output
    if ( $extract->genotypes_output == 0 ) {
        $self->error_message('Executed extract command to create original genotype VCF file, but no genotypes were output. This means they were filtered or ignored because of ambiguous position.');
        return;
    }
    # Check file exists
    if ( not -e $original_genotype_vcf_file ) {
        $self->error_message('Executed extract command to create original genotype VCF file and genotypes were output, but file is gone!');
        return;
    }
    # Check that there are alleles
    my @alleles = grep { $_ ne '--' } keys %{$extract->alleles};
    if ( not @alleles ) {
        $self->error_message('Executed command to create original genotype file, but there are no alleles!');

        return;
    }
    $self->debug_message('Create original genotype VCF file...OK');

    # Original genotype file - has all the info and headers, with positions from this dbsnp
    $self->debug_message('Create original genotype file...');
    my $original_genotype_file = $build->original_genotype_file_path;
    $self->debug_message('Original genotype file: '.$original_genotype_file);
    $extract = Genome::Model::GenotypeMicroarray::Command::Extract->create(
        instrument_data => $instrument_data,
        variation_list_build => $dbsnp_build,
        output => $original_genotype_file.':separator=TAB:fields=chromosome,position,alleles,id,sample_name,log_r_ratio,gc_score,cnv_value,cnv_confidence,allele1,allele2:print_headers=1',
    );
    if ( not $extract ) {
        $self->error_message('Failed to create command to create original genotype file!');
        return;
    }
    $extract->dump_status_messages(1);
    if ( not $extract->execute ) {
        $self->error_message('Failed to execute command to create original genotype file!');
        return;
    }
    # Check that genotypes were output
    if ( $extract->genotypes_output == 0 ) {
        $self->error_message('Executed command to create original genotype file, but no genotypes were output. This means they were filtered or ignored because of ambiguous position.');
        return;
    }
    # Check file exists
    if ( not -e $original_genotype_file ) {
        $self->error_message('Executed command to create original genotype file and genotypes were output, but file is gone!');
        return;
    }
    $self->debug_message('Create original genotype file...OK');

    # Filters for extracting from the above original file
    my @filters = (qw/ gc_score:min=0.7 /); 
    push @filters, 'invalid_iscan_ids' if $reference_sequence_build->version eq '36';

    # Genotype file. No headers, tab sep with chrom, pos and alleles
    $self->debug_message('Create genotype file...');
    my $genotype_file = $build->genotype_file_path;
    $self->debug_message('Genotype file: '.$genotype_file);
    my $extract_genotypes = Genome::Model::GenotypeMicroarray::Command::Extract->create(
        build => $build,
        output => $genotype_file.':separator=TAB:fields=chromosome,position,alleles:print_headers=0',
        filters => \@filters,
    );
    if ( not $extract_genotypes ) {
        $self->error_message('Failed to create command to create genotype file!');
        return;
    }
    $extract_genotypes->dump_status_messages(1);
    if ( not $extract_genotypes->execute ) {
        $self->error_message('Failed to execute command to create genotype file!');
        return;
    }
    if ( not -s $genotype_file ) {
        $self->error_message('Executed command to create genotype file, but file is empty! '.$genotype_file);
        return;
    }
    $self->debug_message('Create genotype file...OK');

    # Nutter made this file name, so we will link to it
    $self->debug_message('Link genotype file to gold2geno file...');
    $self->debug_message('Genotype file: '.$genotype_file);
    my $gold2geno_file = $build->gold2geno_file_path;
    $self->debug_message('Gold2geno file: '.$gold2geno_file);

    # Make a relative symlink if they are in the same directory. I think this
    # will always be the case but since gold2geno_file_path is not locally
    # defined I will check. Relative is better in case build's allocation is
    # moved or archived -> unarchived.
    if (dirname($genotype_file) eq dirname($gold2geno_file)) {
        Genome::Sys->create_symlink(basename($genotype_file), $gold2geno_file);
    } else {
        Genome::Sys->create_symlink($genotype_file, $gold2geno_file);
    }

    if ( not -l $gold2geno_file  or not -s $gold2geno_file ) {
        $self->error_message('Failed to link genotype file to gold2geno file!');
        return;
    }
    $self->debug_message('Link genotype file to gold2geno file...OK');

    # Copy number file. No headers, tab sep with chrom, pos and log r ratio
    $self->debug_message('Create copy number file...');
    my $copy_number_file = $build->copy_number_file_path;
    $self->debug_message('Copy number file: '.$copy_number_file);
    my $extract_copy_number = Genome::Model::GenotypeMicroarray::Command::Extract->create(
        build => $build,
        output => $copy_number_file.':separator=TAB:fields=chromosome,position,log_r_ratio:print_headers=0',
        filters => \@filters,
    );
    if ( not $extract_copy_number ) {
        $self->error_message('Failed to create command to create copy number file!');
        return;
    }
    $extract_copy_number->dump_status_messages(1);
    if ( not $extract_copy_number->execute ) {
        $self->error_message('Failed to execute command to create copy number file!');
        return;
    }
    if ( not -s $copy_number_file ) {
        $self->error_message('Executed command to create copy number file, but file is empty! '.$copy_number_file);
        return;
    }
    $self->debug_message('Create copy number file...OK');

    # TODO bdericks: I'm guessing that second genotype file is supposed to be the replicate. It should be changed
    # to be the actual replicate when we know how to figure it out.
    # abrummet: This is the only place in the tree where this Command is used.  I've stripped out the second input
    # file to fix a bug where it would not read from the "second" file when switching chromosomes and the next position
    # is numerically higher than the last position
    my $snp_array_file = $build->formatted_genotype_file_path;
    $self->debug_message("Create snp array (gold) file: ".$snp_array_file);
    my $gold_snp = Genome::Model::GenotypeMicroarray::Command::CreateGoldSnpFileFromGenotypes->create(
        genotype_file => $genotype_file,
        output_file => $snp_array_file,
        reference_sequence_build => $reference_sequence_build, 
    );
    if ( not $gold_snp ) {
        $self->error_message("Cannot create gold snp tool.");
        return;
    }
    $gold_snp->dump_status_messages(1);
    if ( not $gold_snp->execute ) {
        $self->error_message("Cannot execute gold snp tool");
        return;
    }
    if ( not -s $snp_array_file ) {
        $self->error_message("Executed gold snp tool, but snp array file ($snp_array_file) does not exist");
        return;
    }
    $self->debug_message("Create snp array (gold) file...OK");

    $self->debug_message('Create gold snp bed file...');
    my $snvs_bed = $build->snvs_bed;
    $self->debug_message('Gold snp bed file: '.$snvs_bed);
    my $gold_snp_bed = Genome::Model::GenotypeMicroarray::Command::CreateGoldSnpBed->create(
        input_file => $snp_array_file,
        output_file => $snvs_bed,
        reference => $reference_sequence_build,
    );
    if ( not $gold_snp_bed ) {
        $self->error_message('Failed to create gold snp bed tool!');
        return;
    }
    $gold_snp_bed->dump_status_messages(1);
    unless ($gold_snp_bed->execute) {
        $self->error_message("Could not generate gold snp bed file at $snvs_bed from snp array file $snp_array_file");
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

