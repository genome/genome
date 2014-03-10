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

sub snvs_file {
    my ($self, $version,$format) = @_;
   
    if (not defined($self->version)) {
        $self->error_message("No version set on build?: " . $self->__display_name__);
        # continue for backward compatibility ...should this really work? -sssmith
    }

    # TODO: get a real api for this
    my $name = $self->model->name . "-" . $self->version;
    if (defined $version and $version ne "v1") {
        die "No version of snvs .bed file version $version available for $name";
    }
    
    my $snv_result = $self->snv_result;
    unless ($snv_result) {
        $self->warning_message("No snv result for " . $self->__display_name__);
        return;
    }
    $self->debug_message("Found SNV result: " . $snv_result->__display_name__);

    my $snvs_file_path = join('/', $snv_result->output_dir, "snvs.hq.$format");
    unless (-e $snvs_file_path) {
        $self->error_message("SNVs file not found: $snvs_file_path");
        return;
    }
    $self->debug_message("Found SNV file path: " . $snvs_file_path);

    return $snvs_file_path;
}

1;
