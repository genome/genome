package Genome::InstrumentData::AlignmentResult::Merged::CoverageStats::View::Coverage::Xml;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::AlignmentResult::Merged::CoverageStats::View::Coverage::Xml {
    is => 'UR::Object::View::Default::Xml',
    has => [
        perspective => { default_value => 'coverage' },
    ],
};

#for convenience of overriding in coverage views of other objects
sub coverage_result {
    my $self = shift;

    return $self->subject;
}

sub all_subject_classes {
    my $self = shift;
    my @cl = $self->SUPER::all_subject_classes;

    my $class = 'Genome::InstrumentData::AlignmentResult::Merged::CoverageStats';
    unless(grep( $_ eq $class, @cl )) {
        push @cl, $class;
    }

    return @cl;
}

sub _generate_content {
    my $self = shift;

    my $coverage_result = $self->coverage_result;
    return '' unless $coverage_result;

    my $xml_doc = XML::LibXML->createDocument();
    $self->_xml_doc($xml_doc);

    # the header line is the class followed by the id
    my $object = $xml_doc->createElement('object');
    $xml_doc->setDocumentElement($object);
    my $time = UR::Time->now();
    $object->addChild( $xml_doc->createAttribute("generated-at",$time) );

    $object->addChild( $xml_doc->createAttribute('type', $coverage_result->class) );
    $object->addChild( $xml_doc->createAttribute('id', $coverage_result->id) );

    my $roi_set_name = $object->addChild( $xml_doc->createElement('region_of_interest_set_name') );
    $roi_set_name->addChild( $xml_doc->createTextNode($coverage_result->region_of_interest_set->name) );

    my @instrument_data = $coverage_result->alignment_result->instrument_data;
    my %target_region_set_names;
    for my $i (@instrument_data) {
        $target_region_set_names{$i->target_region_set_name}++;
    }

    my $target_region_set_names = $object->addChild( $xml_doc->createElement('target_region_set_names') );
    for my $target (sort keys %target_region_set_names) {
        $target_region_set_names->addChild( $xml_doc->createTextNode($target) );
    }
    $object->addChild( $self->get_alignment_summary_node );
    $object->addChild( $self->get_coverage_stats_summary_node );
    return $xml_doc->toString(1);
}

sub get_alignment_summary_node {
    my $self = shift;
    my $xml_doc = $self->_xml_doc;

    my $coverage_result = $self->coverage_result;

    my $alignment_summary_hash_ref = $coverage_result->alignment_summary_hash_ref;
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

    my $coverage_result = $self->coverage_result;

    my $coverage_stats_summary_hash_ref = $coverage_result->coverage_stats_summary_hash_ref;
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
