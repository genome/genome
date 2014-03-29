package Genome::Model::Build::GenotypeMicroarray;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::GenotypeMicroarray {
    is => 'Genome::Model::Build',
    has => [
        dbsnp_build_id => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'dbsnp_build', value_class_name => 'Genome::Model::Build::ImportedVariationList' ],
            is_many => 0,
            is_mutable => 1,
            is_optional => 1,
            doc => 'dbsnp build to compare against'
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
    ],
};

sub validate_for_start_methods {
    return (qw/
            validate_inputs_have_values
            inputs_have_compatible_reference
            instrument_data_assigned
            validate_dbsnp_build
    /);
}

sub validate_dbsnp_build {
    my $self = shift;

    my $variation_list_build = $self->dbsnp_build;
    if ( not $variation_list_build ) {
        return UR::Object::Tag->create(
            type => 'error',
            properties => ['dbsnp_build'],
            desc => 'No DB Snp build specified for build!',
        );
    }

    # Need SNVs VCF
    my $snvs_vcf = $variation_list_build->snvs_vcf;
    if ( not defined $snvs_vcf or not -s $snvs_vcf ) {
        return UR::Object::Tag->create(
            type => 'error',
            properties => ['dbsnp_build'],
            desc => 'DB Snp build ('.$variation_list_build->__display_name__.') does not have a SNVS VCF!',
        );
    }

    # Make sure dbsnp has reference
    if ( not $variation_list_build->reference ) {
        return UR::Object::Tag->create(
            type => 'error',
            properties => ['reference_sequence_build'],
            desc => 'DB Snp build does not have a reference sequence build!',
        );
    }

    return;
}

sub perform_post_success_actions {
    my $self = shift;
    return $self->model->request_builds_for_dependent_cron_ref_align;
}

sub create_gold2geno_file_from_genotype_file {
    my $self = shift;

    my $genotype_file = $self->formatted_genotype_file_path;
    my $gold2geno_file = $self->gold2geno_file_path;

    return 1 if (-e $gold2geno_file && $self->validate_gold2geno_file);

    my $genotype_reader = Genome::Sys->open_file_for_reading($genotype_file);
    my $gold2geno_writer = Genome::Sys->open_file_for_writing($gold2geno_file);
    while (my $line = $genotype_reader->getline) {
        my @field = split("\t", $line);
        if ($field[1] ne $field[2]) {
            die $self->error_message("Sample ID differs in Gold SNP file: " . $field[1] . " vs. " . $field[2]);
        }
        $gold2geno_writer->print($field[0] . "\t" . $field[1] . "\t" . $field[3] . $field[4] . "\n");
    }
    if ( -e $gold2geno_file ) {
        $self->validate_gold2geno_file;
    } else {
        die $self->error_message("gold2geno file is missing after conversion.");
    }

    return 1;
}

sub validate_gold2geno_file {
    my $self = shift;

    my $genotype_file = $self->formatted_genotype_file_path;
    my $gold2geno_file = $self->gold2geno_file_path;

    my ($genotype_file_line_count) = qx(wc -l $genotype_file) =~ /^(\d+)/;
    my ($gold2geno_file_line_count) = qx(wc -l $gold2geno_file) =~ /^(\d+)/;

    my $valid_gold2geno_file = ($genotype_file_line_count == $gold2geno_file_line_count);
    if ($valid_gold2geno_file) {
        return 1;
    } else {
        die $self->error_message('gold2geno file exists but line count does not match formatted_genotype_file_path');
    }
}

sub gold2geno_file_path {
    shift->formatted_genotype_file_path . '.gold2geno';
}

sub original_genotype_file_path {
    my $self = shift;
    return $self->data_directory.'/'.$self->model->subject->id.'.original';
}

sub original_genotype_vcf_file_path {
    my $self = shift;
    return $self->original_genotype_file_path.'.vcf';
}

sub genotype_file_path {
    my $self = shift;
    return $self->data_directory.'/'.$self->model->subject->id.'.genotype';
}

sub copy_number_file_path {
    return $_[0]->data_directory.'/'.$_[0]->model->subject->id.'.copynumber';
}

sub formatted_genotype_file_path { # gold snp
    shift->data_directory . '/formatted_genotype_file_path.genotype';
}

sub snvs_bed {
    shift->data_directory . '/gold_snp.v2.bed';
}

sub filtered_snvs_bed {
    shift->data_directory . '/gold_snp.v2.bed';
}

1;

