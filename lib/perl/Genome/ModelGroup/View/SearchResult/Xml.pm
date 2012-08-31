package Genome::ModelGroup::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'id',
                'name',
                {
                    name => 'models',
                    aspects => [ 'id', 'name' ],
                    perspective => 'default',
                    toolkit => 'xml'
                },
                {
                    name => 'convergence_model',
                    aspects => [
                        'id',
                        'name',
                        {
                            name => 'last_complete_build',
                            aspects => [
                                'id', 'data_directory'
                            ],
                            perspective => 'default',
                            toolkit => 'xml',
                            subject_class_name => 'Genome::Model::Build',
                        }
                    ],
                    perspective => 'default',
                    toolkit => 'xml',
                },
            ]
        }
    ]
};

sub _generate_details_result_xml {
    my $self = shift;
    my $subject = $self->subject;
    my $result_node = $self->_result_node;
    
    my $xml_doc = $result_node->ownerDocument;
    
    my $models_node = $result_node->addChild( $xml_doc->createElement('models') );
    
    for my $model ($subject->models) {
        my $model_node = $models_node->addChild( $xml_doc->createElement('model') );
        
        $model_node->addChild( $xml_doc->createAttribute('id', $model->id ));
        $model_node->addChild( $xml_doc->createAttribute('name', $model->name ));
    }
}

1;
