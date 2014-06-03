package Genome::Annotation::Interpreter::FpkmInterpreter;

use strict;
use warnings;
use Genome;

class Genome::Annotation::Interpreter::FpkmInterpreter {
    is => ['Genome::Annotation::Interpreter::Base', 'Genome::Annotation::WithSampleName'],
};

sub name {
    return 'fpkm';
}

sub requires_experts {
    ('fpkm');
}

sub available_fields {
    return qw/
        fpkm
    /;
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;

    my $passed_alt_alleles = shift;

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele} = {
            # FIXME FPKM will have more than one entry depending on the ALT / I wish File::Vcf::Entry supported sample_field_for_alt like it does info fields with info_for_allele -- perhaps add this?
            fpkm => $entry->sample_field($self->sample_index($entry->{header}), "FPKM"),
        };
    }
    return %return_values;
}

1;

