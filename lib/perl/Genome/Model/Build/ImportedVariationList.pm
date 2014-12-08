package Genome::Model::Build::ImportedVariationList;


use strict;
use warnings;

use Data::Dumper;
use Genome;

class Genome::Model::Build::ImportedVariationList {
    is => 'Genome::Model::Build',
    has => [
        version => { 
            via => 'inputs',
            is => 'Text',
            to => 'value_id', 
            where => [ name => 'version', value_class_name => 'UR::Value'], 
            is_mutable => 1 
        },
        reference_id => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'reference', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence' ],
            is_many => 0,
            is_mutable => 1,
            doc => 'reference sequence to align against'
        },
        reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            id_by => 'reference_id',
        },
        source_name => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'source_name', value_class_name => 'UR::Value::Text' ],
            is_many => 0,
            is_mutable => 1,
            is_optional => 1,
            doc => 'The name of the source of the imported variants (e.g., dbsnp, 1kg-wgs)',
        }
    ],
    has_optional_mutable => {
        snv_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            doc => 'The result for snvs to import',
            via => 'inputs',
            to => 'value',
            where => [
                name => 'snv_result',
            ],
        },
        indel_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            doc => 'The result for indels to import',
            via => 'inputs',
            to => 'value',
            where => [
                name => 'indel_result',
            ],
        },
        sv_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            doc => 'The result for svs to import',
            via => 'inputs',
            to => 'value',
            where => [
                name => 'sv_result',
            ],
        },
        cnv_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            doc => 'The result for cnvs to import',
            via => 'inputs',
            to => 'value',
            where => [
                name => 'cnv_result',
            ],
        },
    },
};

sub snvs_bed {
    my ($self, $version) = @_;
    return $self->snvs_file($version, "bed");
}

sub snvs_vcf {
    my ($self, $version) = @_;
    return $self->snvs_file($version, "vcf");
}

sub indels_vcf {
    my ($self, $version) = @_;
    return $self->indels_file($version, "vcf");
}

sub snvs_file {
    my ($self, $version, $format) = @_;
    return $self->file_of_type('snv', $version, $format);
}

sub indels_file {
    my ($self, $version, $format) = @_;
    return $self->file_of_type('indel', $version, $format);
}

sub file_of_type {
    my ($self, $type, $version, $format) = @_;
   
    if (not defined($self->version)) {
        $self->error_message("No version set on build?: " . $self->__display_name__);
        # continue for backward compatibility ...should this really work? -sssmith
    }

    # TODO: get a real api for this
    my $name = $self->model->name . "-" . $self->version;
    if (defined $version and $version ne "v1") {
        die "No version of $type .bed file version $version available for $name";
    }

    my $result_accessor = join("_", $type, "result");
    my $result = $self->$result_accessor;
    unless ($result) {
        $self->warning_message("No $type result for " . $self->__display_name__);
        return;
    }
    $self->debug_message("Found $type result: " . $result->__display_name__);

    my $file_name = $type."s.hq.$format";
    my $file_path = join('/', $result->output_dir, $file_name);
    unless (-e $file_path) {
        $self->error_message("$type file not found: $file_path");
        return;
    }
    $self->debug_message("Found $type file path: " . $file_path);

    return $file_path;
}

1;
