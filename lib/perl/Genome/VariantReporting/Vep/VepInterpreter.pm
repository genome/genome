package Genome::VariantReporting::Vep::VepInterpreter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::VepConsequenceParser;

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

sub process_interpret_entry {
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

1;
