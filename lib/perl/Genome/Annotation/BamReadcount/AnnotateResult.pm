package Genome::Annotation::BamReadcount::AnnotateResult;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::BamReadcount::AnnotateResult {
    is => 'Genome::Annotation::ResultBase',
    has_input => [
        readcount_results => {
            is => 'Genome::Annotation::BamReadcount::RunResult',
            is_many => 1,
        },
    ],
};

sub output_filename {
    return 'with_read_counts.vcf.gz';
}

sub _run {
    my $self = shift;

    Genome::Model::Tools::Vcf::AnnotateWithReadcounts->execute(
        vcf_file => $self->input_result->output_file_path,
        readcount_file_and_sample_idx => $self->readcount_file_and_sample_idxs,
        output_file => File::Spec->join($self->temp_staging_directory, $self->output_filename),
    );

    return;
}

sub readcount_file_and_sample_idxs {
    my $self = shift;
    my $reader = Genome::File::Vcf::Reader->new($self->input_result->output_file_path);
    my $header = $reader->header;

    my @array;
    for my $result ($self->readcount_results) {
        push @array, sprintf("%s:%s", $result->output_file_path,
            $header->index_for_sample_name($result->sample_name));
    }
    return \@array;
}
