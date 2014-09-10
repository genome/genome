package Genome::VariantReporting::Vep::VepInterpreter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::VepConsequenceParser;
use feature "state";

class Genome::VariantReporting::Vep::VepInterpreter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
};

sub name {
    return 'vep';
}

sub requires_annotations {
    return ('vep');
}

sub available_fields {
    return qw/
        transcript_name
        trv_type
        trv_type_category
        amino_acid_change
        default_gene_name
        ensembl_gene_id
        c_position
        gene_name_source
        canonical
        sift
        polyphen
        condel
    /;
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my $vep_parser = new Genome::File::Vcf::VepConsequenceParser($entry->{header});

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        my ($transcript) = $vep_parser->transcripts($entry, $variant_allele);
        $return_values{$variant_allele} = {
            transcript_name   =>$transcript->{'feature'},
            trv_type          => $transcript->{'consequence'},
            trv_type_category => trv_type_category($transcript->{'consequence'}),
            amino_acid_change => $transcript->{'hgvsp'},
            default_gene_name => $transcript->{'symbol'} ,
            ensembl_gene_id   => $transcript->{'gene'},
            gene_name_source => $transcript->{'symbol_source'},
            c_position => $transcript->{'hgvsc'},
            sift => $transcript->{'sift'},
            polyphen => $transcript->{'polyphen'},
            condel => $transcript->{'condel'},
        };
        if ($transcript->{'canonical'} eq "YES") {
            $return_values{$variant_allele}->{'canonical'} = 1,
        }
        else {
            $return_values{$variant_allele}->{'canonical'} = 0,
        }
    }

    return %return_values;
}

sub trv_type_category {
    my $trv_type = shift;

    my $trv_types = Set::Scalar->new(split('&', $trv_type));
    if (is_splice_site($trv_types)) {
        return 'splice_site';
    }
    elsif (is_non_synonymous($trv_types)) {
        return 'non_synonymous';
    }
    else {
        return 'other';
    }
}

sub is_splice_site {
    my $trv_types = shift;

    state $splice_sites = Set::Scalar->new(
        'splice_acceptor_variant',
        'splice_donor_variant'
    );
    return !$splice_sites->intersection($trv_types)->is_null;
}

sub is_non_synonymous {
    my $trv_types = shift;

    state $non_synonymous = Set::Scalar->new(
        'transcript_ablation',
        'stop_gained',
        'frameshift_variant',
        'stop_lost',
        'initiator_codon_variant',
        'transcript_amplification',
        'inframe_insertion',
        'inframe_deletion',
        'missense_variant',
        'incomplete_terminal_codon_variant'
    );
    return !$non_synonymous->intersection($trv_types)->is_null;
}

1;
