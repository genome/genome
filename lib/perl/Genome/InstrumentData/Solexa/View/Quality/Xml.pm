package Genome::InstrumentData::Solexa::View::Quality::Xml;

use strict;
use warnings;
use Genome;
use Data::Dumper;
use XML::LibXML;
use XML::LibXSLT;

class Genome::InstrumentData::Solexa::View::Quality::Xml {
    is => 'UR::Object::View::Default::Xml',
        has => [
            _doc => {
                is_transient => 1,
                doc => 'the XML::LibXML document object used to build the content for this view'
            },
        ]
    };

sub _generate_content {
    my ($self) = @_;

    my $subject = $self->subject;
    return unless $subject;

    my $report_file = get_report("quality.report.xml", $subject);
    my $report_xml = $report_file->content();

    # parse report XML
    my $parser = XML::LibXML->new();
    my $doc = $parser->parse_string( $report_xml );

    $self->_doc($doc);

    # add flow cell ID and other important InstrumentData::Solexa attributes to report node
    my @report_nodes = $doc->getElementsByTagName("report");

    if (scalar(@report_nodes) > 1) { Carp::carp("Found two report nodes when one expected! Ignoring second report."); };

    my $rnode = $report_nodes[0];

    $rnode->addChild( $doc->createAttribute("flow_cell_id", $subject->flow_cell_id) );
    $rnode->addChild( $doc->createAttribute("lane", $subject->lane) );
    $rnode->addChild( $doc->createAttribute("run_name", $subject->run_name) );
    $rnode->addChild( $doc->createAttribute("read_length", $subject->read_length) );
    $rnode->addChild( $doc->createAttribute("sample_id", $subject->sample_id) );
    $rnode->addChild( $doc->createAttribute("sample_name", $subject->sample_name) );

    $rnode->addChild( $doc->createAttribute("filt_error_rate_avg", $subject->filt_error_rate_avg) );

    $rnode->addChild( $doc->createAttribute("gerald_directory", $subject->gerald_directory) );
    $rnode->addChild( $doc->createAttribute("analysis_software_version", $subject->analysis_software_version) );
    $rnode->addChild( $doc->createAttribute("clusters", $subject->clusters) );

    $rnode->addChild( $doc->createAttribute("flow_cell_id", $subject->flow_cell_id) );
    $rnode->addChild( $doc->createAttribute("library_name", $subject->library_name) );

    return $doc->toString(1);

}

sub get_report {
    my ($report_name, $instrument_data) = @_;

    my @fs_ids = GSC::Sequence::ItemFile->get( seq_id => $instrument_data->seq_id );
    my @file = GSC::FileStorage->get( file_storage_id => \@fs_ids, file_name => $report_name  );

    my ($last_file) = sort { $b->file_storage_id <=> $a->file_storage_id } @file;

    return $last_file;

}

1;
