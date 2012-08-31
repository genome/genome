package Genome::DruggableGene::GeneNameGroup::View::Status::Xml;

use strict;
use warnings;
use Genome;
use Data::Dumper;
use XML::LibXML;

class Genome::DruggableGene::GeneNameGroup::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'name',
                {
                    name => 'genes',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'name',
                        'source_db_name',
                        'source_db_version',
                        'source_db_url',
                        'original_data_source_url',
                        {
                            name => 'gene_alt_names',
                            perspective => 'default',
                            toolkit => 'xml',
                            aspects => [
                                'alternate_name',
                                'nomenclature',
                            ],
                        },
                        {
                            name => 'citation',
                            perspective => 'default',
                            toolkit => 'xml',
                            aspects => [
                                'citation',
                            ],
                        },
                    ],
                },
            ],
        },
    ],
};

sub _generate_content {
    my $self = shift;
    my $group = $self->subject;
    my @ids = map{$_->id}$group->genes;
    Genome::DruggableGene::GeneAlternateNameReport->get(gene_id => [@ids]);
    return $self->SUPER::_generate_content(@_);
}

1;
