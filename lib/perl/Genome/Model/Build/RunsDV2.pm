package Genome::Model::Build::RunsDV2;

use strict;
use warnings FATAL => 'all';
use Cwd qw(abs_path);
use Params::Validate qw(validate validate_pos :types);

class Genome::Model::Build::RunsDV2 {
    is_abstract => 1,
};

sub get_detailed_indels_vcf {
    my $self = shift;
    return $self->_get_variant_file("indels.detailed.vcf.gz");
}

sub get_detailed_snvs_vcf {
    my $self = shift;
    return $self->_get_variant_file("snvs.detailed.vcf.gz");
}

sub get_detailed_pass_indels_vcf {
    my $self = shift;
    return $self->_get_variant_file("indels.detailed.pass.vcf.gz");
}

sub get_detailed_pass_snvs_vcf {
    my $self = shift;
    return $self->_get_variant_file("snvs.detailed.pass.vcf.gz");
}

sub get_indels_vcf {
    my $self = shift;
    return $self->_get_variant_file("indels.vcf.gz");
}

sub get_snvs_vcf {
    my $self = shift;
    return $self->_get_variant_file("snvs.vcf.gz");
}

sub _get_variant_file {
    my ($self, $filename) = @_;
    return File::Spec->join($self->variants_directory, $filename);
}

sub variants_directory {
    my $self = shift;

    my $expected_directory = File::Spec->join($self->data_directory, 'variants');

    if (-d $expected_directory) {
        return $expected_directory;
    } else {
        #for compatibility with previously existing builds
        for my $dir_name (qw(snp_related_metrics  sam_snp_related_metrics  maq_snp_related_metrics  var-scan_snp_related_metrics)) {
            my $dir = File::Spec->join($self->data_directory, $dir_name);
            return $dir if -d $dir;
        }
    }
    return $expected_directory;
}

sub get_detailed_vcf_result {
    my ($self, $type) = validate_pos(@_, 1, 1);
    my $accessor = sprintf('get_detailed_%s_vcf_result', $type);
    return $self->$accessor;
}

sub get_detailed_indels_vcf_result {
    my $self = shift;
    return _get_result_for_file($self->get_detailed_indels_vcf());
}

sub get_detailed_snvs_vcf_result {
    my $self = shift;
    return _get_result_for_file($self->get_detailed_snvs_vcf());
}

sub _get_result_for_file {
    my $path = shift;
    my $allocation = Genome::Disk::Allocation->get_allocation_for_path(abs_path($path));
    unless (defined $allocation) {
        die sprintf("The allocation for path %s doesn't exist", $path);
    }

    my $candidate = $allocation->owner;
    my $expected_isa = 'Genome::Model::Tools::DetectVariants2::Result::Vcf';
    if (defined $candidate and $candidate->isa($expected_isa)) {
        return $candidate;
    } else {
        die sprintf("The owner (%s) of allocation (%s) for path (%s) isn't a subclass of (%s) as expected.",
            ref($candidate), $allocation->id, $path, $expected_isa);
    }
}

sub _dv2_result_subclass_names {
    my $class = shift;

    return qw(
        Genome::Model::Tools::DetectVariants2::Classify::Loh
        Genome::Model::Tools::DetectVariants2::Classify::PreviouslyDiscovered
        Genome::Model::Tools::DetectVariants2::Classify::Tier
        Genome::Model::Tools::DetectVariants2::Result
        Genome::Model::Tools::DetectVariants2::Result::Combine::IntersectSnv
        Genome::Model::Tools::DetectVariants2::Result::Combine::IntersectIndel
        Genome::Model::Tools::DetectVariants2::Result::Combine::LqUnion
        Genome::Model::Tools::DetectVariants2::Result::Combine::UnionCnv
        Genome::Model::Tools::DetectVariants2::Result::Combine::UnionSv
        Genome::Model::Tools::DetectVariants2::Result::Combine::UnionIndel
        Genome::Model::Tools::DetectVariants2::Result::Combine::UnionuniqueIndel
        Genome::Model::Tools::DetectVariants2::Result::Combine::UnionuniqueSnv
        Genome::Model::Tools::DetectVariants2::Result::Filter
        Genome::Model::Tools::DetectVariants2::Result::Vcf::Combine
        Genome::Model::Tools::DetectVariants2::Result::Vcf::Detector
        Genome::Model::Tools::DetectVariants2::Result::Vcf::Filter
    );
}

1;
