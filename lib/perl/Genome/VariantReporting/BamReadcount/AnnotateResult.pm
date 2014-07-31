package Genome::VariantReporting::BamReadcount::AnnotateResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::File::Vcf::Reader;

class Genome::VariantReporting::BamReadcount::AnnotateResult {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Result',
    has_input => [
        readcount_results => {
            is => 'Genome::VariantReporting::BamReadcount::RunResult',
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
        vcf_file => $self->input_vcf,
        readcount_file_and_sample_name => $self->readcount_file_and_sample_names,
        output_file => File::Spec->join($self->temp_staging_directory, $self->output_filename),
    );

    return;
}

sub readcount_file_and_sample_names {
    my $self = shift;
    my $reader = Genome::File::Vcf::Reader->new($self->input_vcf);
    my $header = $reader->header;

    my @array;
    for my $result ($self->readcount_results) {
        push @array, sprintf("%s:%s", $result->output_file_path,
            $result->sample_name);
    }
    return \@array;
}
