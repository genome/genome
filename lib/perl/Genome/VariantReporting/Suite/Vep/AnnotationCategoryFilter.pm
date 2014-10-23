package Genome::VariantReporting::Suite::Vep::AnnotationCategoryFilter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::VepConsequenceParser;
use Genome::VariantReporting::Suite::Vep::AnnotationCategory;

class Genome::VariantReporting::Suite::Vep::AnnotationCategoryFilter {
    is  => ['Genome::VariantReporting::Framework::Component::Filter'],
    has => {
        category_list => {
            is => 'Text',
            is_many => 1,
        },
    },
};


sub name {
    return 'annotation-category';
}

sub requires_annotations {
    return ('vep');
}


sub filter_entry {
    my ($self, $entry) = @_;
    my %return_values;

    my $vep_parser = Genome::File::Vcf::VepConsequenceParser->new($entry->{header});

    ALLELE: for my $alt_allele (@{$entry->{alternate_alleles}}) {
        my ($transcript) = $vep_parser->transcripts($entry, $alt_allele);
        my $consequence  = $transcript->{consequence};
        my $value        = 0;

        if (defined $consequence) {
            my @types = split /\&/, lc($consequence);

            LIST: for my $category ($self->category_list) {
                my $method = 'is_'.$category;
                my $category_class = 'Genome::VariantReporting::Suite::Vep::AnnotationCategory';

                unless ($category_class->can($method)) {
                    die $self->error_message('No method (%s) defined for (%s) in %s', $method, $category, $category_class);
                }

                if ($category_class->$method(@types)) {
                    $value = 1;
                    last LIST;
                }
            }
        }
        $return_values{$alt_allele} = $value;
    }

    return %return_values;
}


sub vcf_id {
    my $self = shift;
    my @list = map{s/\_//g}$self->category_list;
    return 'ANNOT_' . join '_', @list;
}


sub vcf_description {
    my $self = shift;
    return 'Variant hits annotation categories: ' . join ',', $self->category_list;
}

1;




