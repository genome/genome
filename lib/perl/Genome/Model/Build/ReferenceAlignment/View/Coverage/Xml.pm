package Genome::Model::Build::ReferenceAlignment::View::Coverage::Xml;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::ReferenceAlignment::View::Coverage::Xml {
    is => 'Genome::InstrumentData::AlignmentResult::Merged::CoverageStats::View::Coverage::Xml',
};

sub coverage_result {
    my $self = shift;

    my $build = $self->subject;
    return unless $build;

    return $build->reference_coverage_result;
}

sub _generate_content {
    my $self = shift;

    if($self->coverage_result) {
        #new view for result based builds
        return $self->SUPER::_generate_content;
    }

    #old view for non-result based builds

    my $subject = $self->subject();
    return '' unless $subject;

    my $model = $subject->model;

    my $xml_doc = XML::LibXML->createDocument();
    $self->_xml_doc($xml_doc);

    # the header line is the class followed by the id
    my $object = $xml_doc->createElement('object');
    $xml_doc->setDocumentElement($object);
    my $time = UR::Time->now();
    $object->addChild( $xml_doc->createAttribute("generated-at",$time) );

    $object->addChild( $xml_doc->createAttribute('type', $self->subject_class_name) );
    $object->addChild( $xml_doc->createAttribute('id', $subject->id) );

    my $roi_set_name = $object->addChild( $xml_doc->createElement('region_of_interest_set_name') );
    $roi_set_name->addChild( $xml_doc->createTextNode($model->region_of_interest_set_name) );

    my $target_region_set_names = $object->addChild( $xml_doc->createElement('target_region_set_names') );
    my @inputs = Genome::Model::Input->get(model_id => $model->id, name => 'target_region_set_name');
    for my $input (@inputs) {
        $target_region_set_names->addChild( $xml_doc->createTextNode($input->value_id) );
    }
    $object->addChild( $self->get_alignment_summary_node );
    $object->addChild( $self->get_coverage_stats_summary_node );
    return $xml_doc->toString(1);
}

sub get_alignment_summary_node {
    my $self = shift;
    my $xml_doc = $self->_xml_doc;
    my $build = $self->subject;
    my $alignment_summary_hash_ref = $build->alignment_summary_hash_ref;
    my $alignment_summary_node = $xml_doc->createElement('alignment-summary');
    for my $wingspan (keys %{$alignment_summary_hash_ref}) {
        my $wingspan_node = $alignment_summary_node->addChild( $xml_doc->createElement('wingspan') );
        $wingspan_node->addChild( $xml_doc->createAttribute('value',$wingspan) );
        for my $key (keys %{$alignment_summary_hash_ref->{$wingspan}}) {
            my $key_node = $wingspan_node->addChild( $xml_doc->createElement($key) );
            $key_node->addChild( $xml_doc->createTextNode( $alignment_summary_hash_ref->{$wingspan}->{$key} ) );
        }
    }
    return $alignment_summary_node;
}

sub get_coverage_stats_summary_node {
    my $self = shift;
    my $xml_doc = $self->_xml_doc;
    my $build = $self->subject;
    my $coverage_stats_summary_hash_ref = $build->coverage_stats_summary_hash_ref;
    my $coverage_stats_summary_node = $xml_doc->createElement('coverage-stats-summary');
    for my $wingspan (keys %{$coverage_stats_summary_hash_ref}) {
        my $wingspan_node = $coverage_stats_summary_node->addChild( $xml_doc->createElement('wingspan') );
        $wingspan_node->addChild( $xml_doc->createAttribute('value',$wingspan) );
        for my $min_depth (keys %{$coverage_stats_summary_hash_ref->{$wingspan}}) {
            my $min_depth_node = $wingspan_node->addChild( $xml_doc->createElement('minimum_depth') );
            $min_depth_node->addChild( $xml_doc->createAttribute('value',$min_depth) );
            for my $key (keys %{$coverage_stats_summary_hash_ref->{$wingspan}->{$min_depth}}) {
                my $key_node = $min_depth_node->addChild( $xml_doc->createElement($key) );
                $key_node->addChild( $xml_doc->createTextNode( $coverage_stats_summary_hash_ref->{$wingspan}->{$min_depth}->{$key} ) );
            }
        }
    }
    return $coverage_stats_summary_node;
}
