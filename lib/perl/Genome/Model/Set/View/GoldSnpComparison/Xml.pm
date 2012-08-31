package Genome::Model::Set::View::GoldSnpComparison::Xml;

use strict;
use warnings;

use Genome;

class Genome::Model::Set::View::GoldSnpComparison::Xml {
    is => 'UR::Object::View::Default::Xml',
    has_constant => [
        perspective => {
            value => 'gold-snp-comparison',
        },
    ]
};

sub _generate_content {
    my $self = shift;
    my $subject = $self->subject;
    
    return unless $subject;
    
    my @models = $self->_subject_models;

    my $doc = XML::LibXML->createDocument();
    $self->_xml_doc($doc);
    my $report_node = $doc->createElement("report");
    $doc->setDocumentElement($report_node);
    $report_node->addChild($self->_get_meta_node());

    my $dataset_node = $doc->createElement("dataset");

    for my $model (@models) {
        $dataset_node->addChild($self->_get_model_node($model));
    }
    $report_node->addChild($dataset_node);
    
    return $self->_xml_doc->toString(1);
}

sub _subject_models {
    my $self = shift;
    my @models = $self->subject->members;

}


sub _get_model_node {
    my $self = shift;
    my $model = shift;

    my $model_node = $self->anode("model","id",$model->id,"name",$model->name,"processing-profile",$model->processing_profile_name,"username",$model->user_name);

    my @builds = $model->builds;
    for (@builds) {
        my $metrics = $self->_metrics_for_build($_);
        next if ($metrics->{gigabases} == 0);

        my $build_node = $self->anode("build","id",$_->id);
        $build_node->addChild($self->tnode('lanes',$metrics->{lanes}));
        $build_node->addChild($self->tnode('gigabases',$metrics->{gigabases}));
        $build_node->addChild($self->tnode('goldsnp-concordance-filtered',$metrics->{goldsnp_concordance_filtered}));

        $model_node->addChild($build_node);
    }
   
    return $model_node; 
}

sub _get_meta_node {
    my $self =  shift;

    my @models = $self->_subject_models;

    my $meta_node = $self->_xml_doc->createElement("report-meta");

    $meta_node->addChild($self->tnode('name','Model Builds Report'));
    $meta_node->addChild($self->tnode('description','Shows the gigabases, filtered gold SNP concordance, and flow cells for all the builds in a model'));
    $meta_node->addChild($self->tnode('date',UR::Time->now()));
    $meta_node->addChild($self->tnode('generator',ref($self)));
    
    my $params_node = $self->_xml_doc->createElement("generator-params");
    for my $model_node (map {$self->tnode('model-id', $_->id)} @models) {
        $params_node->addChild($model_node);
    }
    $meta_node->addChild($params_node);

    return $meta_node;

}

sub _metrics_for_build {
    my $self = shift;
    my $build = shift;
    
    my $het_metric_name='gold-heterozygous-snp match heterozygous-1-allele-variant filtered';
    my $kb_metric_name='instrument data total kb';

    my $het_value = ($build->get_metric($het_metric_name) || 0);
    my $kb_value = ($build->get_metric($kb_metric_name) || 0);
    
    my $gb_value = $kb_value/1000000;

    my @instrument_data = $build->instrument_data;
    return {"lanes"=>scalar @instrument_data, "goldsnp_concordance_filtered" => $het_value, "gigabases" => $gb_value};
}

sub create_node_with_attribute {
    
    my $self = shift;
    my $node_name = shift;

    my $doc = $self->_xml_doc;
    
    my $node = $doc->createElement($node_name);
    while (my $attr_name = shift ) {
        my $attr_value = shift;
        $node->addChild($doc->createAttribute($attr_name,$attr_value));
    }
    return $node;
} 

#helper methods.  just pass through to the more descriptive names 
#anode = attribute node
sub anode {
    my $self = shift;
    return $self->create_node_with_attribute(@_);
}
 
#tnode = text node
sub tnode {
    my $self = shift; 
    return $self->create_node_with_text(@_);
}
 
sub create_node_with_text {
    my $self = shift;
    my $node_name = shift;
    my $node_value = shift;

    my $doc = $self->_xml_doc;
    
    my $node = $doc->createElement($node_name);
    if ( defined($node_value) ) {
        $node->addChild($doc->createTextNode($node_value));
    } 
    return $node;

} 

