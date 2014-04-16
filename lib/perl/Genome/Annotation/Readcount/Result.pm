package Genome::Annotation::Readcount::Result;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Readcount::Result {
    is => 'Genome::Annotation::Detail::Result',

    has_input => [
        input_vcf_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Vcf',
        },
        readcount_results => {
            is => 'Genome::Annotation::RunBamReadcount::Result',
            is_many => 1,
        },
    ],
    has_param => [
        variant_type => {
            is => 'Text',
            valid_values => ['snvs', 'indels'],
        },
    ],
};

sub vcf_filename {
    return 'with_read_counts.vcf.gz';
}

sub vcf_file {
    my $self = shift;

    return File::Spec->join($self->output_dir, $self->vcf_filename);
}

sub _run {
    my $self = shift;

    Genome::Model::Tools::Vcf::AnnotateWithReadcounts->execute(
        vcf_file => $self->input_vcf_file,
        readcount_file_and_sample_idx => $self->readcount_file_and_sample_idxs,
        output_file => File::Spec->join($self->temp_staging_directory, $self->vcf_filename),
    );

    return;
}

sub input_vcf_file {
    my $self = shift;

    return $self->input_vcf_result->get_vcf($self->variant_type),
}

sub readcount_file_and_sample_idxs {
    my $self = shift;
    my $reader = Genome::File::Vcf::Reader->new($self->input_vcf_file);
    my $header = $reader->header;

    my @array;
    for my $result ($self->readcount_results) {
        push @array, sprintf("%s:%s", $result->readcount_file,
            $header->index_for_sample_name($result->sample_name));
    }
    return \@array;
}
