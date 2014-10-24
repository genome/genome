package Genome::VariantReporting::Suite::Vep::VepInterpreter;

use strict;
use warnings;

use Genome;
use Genome::File::Vcf::VepConsequenceParser;
use Genome::VariantReporting::Suite::Vep::AnnotationCategory;

class Genome::VariantReporting::Suite::Vep::VepInterpreter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
};

sub name {
    return 'vep';
}

sub requires_annotations {
    return ('vep');
}

sub field_descriptions {
    return (
        transcript_name   => 'Ensembl stable ID of feature',
        trv_type          => 'Consequence type of this variation',
        trv_type_category => 'The category for the trv type: splice_site, non_synonymous, other',
        amino_acid_change => 'The HGVS protein sequence name',
        default_gene_name => 'The gene symbol',
        ensembl_gene_id   => 'Ensembl stable ID of affected gene',
        c_position        => 'The HGVS coding sequence name',
        gene_name_source  => 'The source of the gene symbol',
        canonical         => 'Is this transcript canonical',
        sift              => 'The SIFT prediction and/or score, with both given as prediction(score)',
        polyphen          => 'The PolyPhen prediction and/or score',
        condel            => 'The Consensus Deleteriousness (Condel) score',
    )
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
            transcript_name   => $transcript->{'feature'},
            trv_type          => $transcript->{'consequence'},
            trv_type_category => trv_type_category($transcript->{'consequence'}),
            amino_acid_change => $transcript->{'hgvsp'},
            default_gene_name => $transcript->{'symbol'} ,
            ensembl_gene_id   => $transcript->{'gene'},
            gene_name_source  => $transcript->{'symbol_source'},
            c_position        => $transcript->{'hgvsc'},
            sift              => $transcript->{'sift'},
            polyphen          => $transcript->{'polyphen'},
            condel            => $transcript->{'condel'},
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
    my $type_info = shift;
    my @types     = split('&', $type_info);
    my $category  = 'Genome::VariantReporting::Suite::Vep::AnnotationCategory';
    my $trv_type  = 'other';

    for my $category_type qw(splice_site non_synonymous) {
        if ($category->is_category($category_type, @types)) {
            $trv_type = $category_type;
            last;
        }
    }
    return $trv_type;
}


1;
