package Genome::Annotation::Readcount::Result;

use strict;
use warnings FATAL => 'all';
use Genome;
use Sys::Hostname;

class Genome::Annotation::Readcount::Result {
    is => 'Genome::SoftwareResult::Stageable',
    has_input => [
        aligned_bam_result => {
            is => 'Genome::InstrumentData::AlignedBamResult',
        },
        vcf_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Vcf',
        },
    ],
    has_param => [
        use_version => { is  => 'Version', },
        variant_type => { is => 'Text', },
        minimum_mapping_quality => { is => 'Integer', },
        minimum_base_quality => { is => 'Integer', },
        max_count => { is  => 'Integer', },
        per_library => { is  => 'Bool', },
        insertion_centric => { is  => 'Bool', },
    ],
};

sub readcount_filename {
    return 'bam-readcount-output.tsv';
}

sub readcount_file {
    my $self = shift;

    return File::Spec->join($self->output_dir, $self->readcount_filename);
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->_prepare_staging_directory;
    $self->_run;

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;
}

sub _run {
    my $self = shift;

    my $region_list = Genome::Sys->create_temp_file_path();
    Genome::Model::Tools::Bed::Convert::VcfToBed->execute(
        remove_filtered_calls => 0,
        source => $self->vcf_result->get_vcf($self->variant_type),
        output => $region_list,
        one_based => 1,
    );

    # Sam::Readcount doesn't accept variant_type
    my %params = $self->param_hash;
    delete $params{variant_type};

    Genome::Model::Tools::Sam::Readcount->execute(
        bam_file => $self->aligned_bam_result->bam_file,
        reference_fasta => $self->aligned_bam_result->reference_fasta,
        output_file => File::Spec->join($self->temp_staging_directory, $self->readcount_filename),
        region_list => $region_list,
        %params,
    );

    return;
}

sub param_names {
    my $self = shift;

    my @properties = $self->__meta__->properties(class_name => __PACKAGE__);
    return map {$_->property_name} grep {$_->is_param} @properties;
}

sub param_hash {
    my $self = shift;

    my %hash;
    for my $param_name ($self->param_names) {
        $hash{$param_name} = $self->$param_name;
    }
    return %hash;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    return File::Spec->join('/', 'model_data', 'software-result', $self->id);
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}
