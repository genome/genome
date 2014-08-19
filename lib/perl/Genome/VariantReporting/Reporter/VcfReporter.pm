package Genome::VariantReporting::Reporter::VcfReporter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::Writer;
use List::AllUtils qw(any);

class Genome::VariantReporting::Reporter::VcfReporter {
    is => 'Genome::VariantReporting::Framework::Component::Reporter::SingleFile',
    has_transient_optional => [
        vcf_file => {
            is => 'Genome::File::Vcf',
        },
        header => {
            is => 'Genome::Vcf::Header',
        },
        header_is_written => {
            is => 'Bool',
        },
    ],
};

sub name {
    return 'vcf';
}

sub requires_interpreters {
    return qw(vcf-entry min-coverage);
}

sub report {
    my $self = shift;
    my $interpretations = shift;

    my %vcf_entry_interpretations = %{$interpretations->{'vcf-entry'}};
    # All the vcf-entry interpretations are duplicates due to the caveat of
    # having to pass interpretations for each passed alt allele
    # We just grab the first one
    my $entry = (values %vcf_entry_interpretations)[0]->{'vcf_entry'};

    unless ($self->header_is_written) {
        my $header = $entry->{'header'};
        $self->add_headers_for_soft_filters($header);
        $self->add_header_for_main_filter($header);
        $self->header($header);
        $self->print_vcf_header();
    }

    my @passed_alleles = keys %vcf_entry_interpretations;
    #get interpretations of filter interpreters on a per allele basis
    #pass allele if all alt alleles are pass
    #add filter status to the info field

    #print the entry
    $self->vcf_file->write($entry);

}

sub print_vcf_header {
    my $self = shift;

    $self->vcf_file(Genome::File::Vcf::Writer->fhopen($self->_output_fh, $self->file_name, $self->header));
    $self->header_is_written(1);
}

sub add_headers_for_soft_filters {
    my $self = shift;
    my $header = shift;

    for my $filter_interpreter ($self->filter_interpreters) {
        $header->add_filter_str(sprintf(
            "<ID=%s,Description=\"%s\">",
            $filter_interpreter->vcf_id,
            $filter_interpreter->vcf_description,
        ));
    }
}

sub add_header_for_main_filter {
    my $self = shift;
    my $header = shift;

    my $filters = join(", ", map { $_->vcf_id } $self->filter_interpreters);

    $header->add_info_str(sprintf(
        "<ID=%s,Number=A,Type=Flag,Description=\"%s\">",
        'ALLFILTERSPASS',
        'Did the alternate allele pass all the soft filters: ' . $filters
    ));
}

sub filter_interpreters {
    my $self = shift;
    return grep { $_->is_filter_interpreter } values %{$self->interpreters};
}

1;
