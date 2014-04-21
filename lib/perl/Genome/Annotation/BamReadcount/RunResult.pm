package Genome::Annotation::BamReadcount::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use Sys::Hostname;

class Genome::Annotation::BamReadcount::RunResult {
    is => 'Genome::Annotation::ResultBase',
    has_input => [
        aligned_bam_result => {
            is => 'Genome::InstrumentData::AlignedBamResult',
        },
    ],
    has_param => [
        version => { is  => 'Version', },
        minimum_mapping_quality => { is => 'Integer', },
        minimum_base_quality => { is => 'Integer', },
        max_count => { is  => 'Integer', },
        per_library => { is  => 'Bool', },
        insertion_centric => { is  => 'Bool', },
    ],
};

sub output_filename {
    return 'bam-readcount-output.tsv';
}

sub _run {
    my $self = shift;

    my $region_list = Genome::Sys->create_temp_file_path();
    Genome::Model::Tools::Bed::Convert::VcfToBed->execute(
        remove_filtered_calls => 0,
        source => $self->input_result->output_file_path,
        output => $region_list,
        one_based => 1,
    );

    # Sam::Readcount doesn't accept variant_type or test_name
    my %params = $self->param_hash;
    delete $params{variant_type};
    delete $params{test_name};

    Genome::Model::Tools::Sam::Readcount->execute(
        bam_file => $self->aligned_bam_result->bam_file,
        reference_fasta => $self->aligned_bam_result->reference_fasta,
        output_file => File::Spec->join($self->temp_staging_directory, $self->output_filename),
        region_list => $region_list,
        use_version => delete $params{version},
        %params,
    );

    return;
}

sub sample_name {
    my $self = shift;

    return $self->aligned_bam_result->sample_name;
}
