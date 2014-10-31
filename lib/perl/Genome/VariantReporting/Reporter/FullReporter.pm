package Genome::VariantReporting::Reporter::FullReporter;

use strict;
use warnings;
use Genome;
use List::Util qw( min );
use Genome::VariantReporting::Suite::BamReadcount::VafInterpreterHelpers qw(single_vaf_headers per_library_vaf_headers);

class Genome::VariantReporting::Reporter::FullReporter {
    is => [ 'Genome::VariantReporting::Reporter::WithHeader', 'Genome::VariantReporting::Framework::Component::WithManySampleNames'],
    has => [
    ],
};

sub name {
    return 'full';
}

sub requires_interpreters {
    return qw(position vep info-tags variant-type min-coverage min-coverage-observed max-vaf-observed variant-callers many-samples-vaf rsid caf);
}

sub headers {
    my $self = shift;
    my @headers = qw/
        chromosome_name
        start
        stop
        reference
        variant
        variant_type
        transcript_name
        trv_type
        trv_type_category
        amino_acid_change
        default_gene_name
        ensembl_gene_id
        inSegDup
        AML_RMG
        rsid
        caf
        max_alt_af
        onTarget
        MeetsMinDepthCutoff
    /;

    push @headers, single_vaf_headers([$self->sample_names]);

    push @headers, qw/
        min_coverage_observed
        max_normal_vaf_observed
        max_tumor_vaf_observed
        variant_callers
        variant_caller_count
    /;

    push @headers, per_library_vaf_headers([$self->sample_names]);

    return @headers;
}

sub _header_to_info_tag_conversion {
    return {
        AML_RMG => 'AML_RMG',
        inSegDup => 'SEG_DUP',
        onTarget => 'ON_TARGET',
    }
}

sub report {
    my $self = shift;
    my $interpretations = shift;

    my %fields = $self->available_fields_dict();
    for my $allele (keys %{$interpretations->{($self->requires_interpreters)[0]}}) {
        for my $header ($self->headers()) {
            my $interpreter = $fields{$header}->{interpreter};
            my $field = $fields{$header}->{field};

            if ($header eq 'inSegDup' || $header eq 'onTarget' || $header eq 'AML_RMG') {
                my $info_tags = $interpretations->{$interpreter}->{$allele}->{$field};
                $self->_print_info_tag(_header_to_info_tag_conversion()->{$header}, $info_tags);
            }
            elsif ($header eq 'variant_callers') {
                my @variant_callers = @{$interpretations->{$interpreter}->{$allele}->{$field}};
                $self->_output_fh->print($self->_format(join(", ", @variant_callers)) . "\t");
            }
            # If we don't have an interpreter that provides this field, handle it cleanly if the field is known unavailable
            elsif ($self->header_is_unavailable($header)) {
                $self->_output_fh->print( $self->_format() . "\t");
            }
            elsif ($interpreter) {
                $self->_output_fh->print($self->_format($interpretations->{$interpreter}->{$allele}->{$field}) . "\t");
            }
            else {
                # We use $header here because $field will be undefined due to it not being in an interpreter
                die $self->error_message("Field (%s) is not available from any of the interpreters provided", $header);
            }
        }
        $self->_output_fh->print("\n");
    }
}
sub available_fields_dict {
    my $self = shift;

    my %available_fields = $self->SUPER::available_fields_dict();
    for my $info_tag_field (qw/inSegDup onTarget AML_RMG/) {
        $available_fields{$info_tag_field} = {
            interpreter => 'info-tags',
            field => 'info_tags',
        };
    }

    $available_fields{MeetsMinDepthCutoff} = {
        interpreter => 'min-coverage',
        field => 'filter_status',
    };

    return %available_fields;
}

sub _print_info_tag {
    my ($self, $info_tag, $info_tags_string) = @_;

    if ($info_tags_string =~ /$info_tag/) {
        $self->_output_fh->print($self->_format(1) . "\t");
    }
    else {
        $self->_output_fh->print($self->_format(0) . "\t");
    }
}
1;
