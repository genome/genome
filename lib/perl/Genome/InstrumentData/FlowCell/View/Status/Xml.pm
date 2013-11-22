package Genome::InstrumentData::FlowCell::View::Status::Xml;

use strict;
use warnings;
use Genome;
use Data::Dumper;
use XML::LibXML;
use XML::LibXSLT;
use Switch;

class Genome::InstrumentData::FlowCell::View::Status::Xml {
    is => 'UR::Object::View::Default::Xml',
        has => [
            _doc    => {
                is_transient => 1,
                doc => 'the XML::LibXML document object used to build the content for this view'
            },
        ],
            has_optional => [
                section => {
                    is => 'String',
                    doc => "NOT IMPLEMENTED YET.  The sub-section of the document to return.  Options are 'all', 'lanes', etc.",
                }
            ],
        };

# this is expected to return an XML string
# it has a "subject" property which is the flowcell we're viewing

sub _generate_content {
    my $self = shift;

    #create the XML doc and add it to the object
    my $doc = XML::LibXML->createDocument();
    $self->_doc($doc);

    my $subject = $self->subject;
    return unless $subject;

    my $flowcell_node = $doc->createElement('flow-cell');
    $flowcell_node->addChild( $doc->createAttribute('id', $subject->flow_cell_id) );

    my $production_node = $flowcell_node->addChild( $doc->createElement('production') );
    $production_node->addChild( $doc->createAttribute('date-started', $subject->production_started) );
    $production_node->addChild( $doc->createAttribute('run-name', $subject->run_name) );
    $production_node->addChild( $doc->createAttribute('run-type', $subject->run_type) );
    $production_node->addChild( $doc->createAttribute('group-name', $subject->group_name) );
    $production_node->addChild( $doc->createAttribute('machine-name', $subject->machine_name) );
    # $production_node->addChild( $doc->createAttribute('team-name', $subject->team_name) );

    if ($subject->lane_info) {
        for my $lane ($subject->lane_info) {
            my $instrument_data_node = $flowcell_node->addChild( $doc->createElement('instrument-data') );
            $instrument_data_node->addChild( $doc->createAttribute('id', ${$lane}{id}) );
            $instrument_data_node->addChild( $doc->createAttribute('lane', ${$lane}{lane}) );
            $instrument_data_node->addChild( $doc->createAttribute('subset_name', ${$lane}{subset_name}) );
            my $gerald_directory = ${$lane}{gerald_directory};
            $instrument_data_node->addChild( $doc->createAttribute('gerald-directory', $gerald_directory) )
                if ($gerald_directory and -e $gerald_directory);;

            for my $url (@{ $lane->{lane_reports} }) {
                # determine report name from report file
                my $name;

                switch($url) {
                    # quality report
                    case m/quality\.html/ { $name = "quality" };

                    # gc_bias report
                    case m/gc-bias-chart\.pdf/ {
                        $name = "gc bias"
                    };

                    #fastqc report
                    case m/fastqc_report/ {
                        $url =~ m/fastqc\/(.*)_sequence_fastqc/;
                        $name = "fast qc (" . $1 . ")";
                    };

                }

                if ($url =~ m/gscmnt/)  { $url = $ENV{GENOME_SYS_SERVICES_FILES_URL} . $url; }

                my $report_node = $instrument_data_node->addChild( $doc->createElement('report'));
                $report_node->addChild( $doc->createAttribute('name', $name) );
                $report_node->addChild( $doc->createAttribute('url', $url) );
            }
        }
    }

    if ($subject->illumina_index) {
        my ($index_sequences, $indexes) = $subject->illumina_index;
        # $index_sequences contains the indexes that are expected for each lane

        my $il_index_node = $flowcell_node->addChild( $doc->createElement('illumina-lane-index') );
        my $kb_report = $il_index_node->addChild( $doc->createElement('report') );
        $kb_report->addChild( $doc->createAttribute('name', 'kilobases_read' ) );

        my $total_read;

        foreach my $lane (keys %{$indexes}) {
            foreach my $index_sequence (keys %{$indexes->{$lane}}) {
                $total_read += $indexes->{$lane}->{$index_sequence}->{'fwd_kilobases_read'};
            }
        }

        foreach my $lane (sort keys %{$indexes}) {
            my $lane_node = $kb_report->addChild( $doc->createElement('lane') );
            $lane_node->addChild( $doc->createAttribute('number', $lane) );

            for my $sequence (sort keys %{$index_sequences}) {
                my $index_node = $lane_node->addChild( $doc->createElement('index') );
                my $sequence_node = $index_node->addChild( $doc->createElement('sequence') );
                $sequence_node->addChild( $doc->createTextNode( $sequence ) );
                my $percent_node = $index_node->addChild( $doc->createElement('percent') );
                my $count_node = $index_node->addChild( $doc->createElement('count') );

                # if we find the expected index sequence, add its info to the index node
                if (defined($indexes->{$lane}->{$sequence})) {
                    my $fwd_bases = $indexes->{$lane}->{$sequence}->fwd_kilobases_read;
                    if (defined($fwd_bases)) {
                        my $percentage;
                        if ($fwd_bases == 0) {
                            $percentage = "0";
                        }
                        elsif ($fwd_bases != 0 and $total_read == 0) {
                            $percentage = "ERROR";
                        }
                        else {
                            $percentage = sprintf("%.2f", ($fwd_bases / $total_read) * 100);
                        }
                        $percent_node->addChild($doc->createTextNode($percentage));
                    } else {
                        $percent_node->addChild( $doc->createTextNode( "0" ) );
                    }
                    $count_node->addChild( $doc->createTextNode($fwd_bases) );
                } else {
                    # no data found for this index in this lane,
                    # so we insert 0 nodes for percent and count
                    $percent_node->addChild( $doc->createTextNode( "0" ) );
                    $count_node->addChild( $doc->createTextNode( "0" ) );
                }
            }
        }
    }

    #set the build status node to be the root
    $doc->setDocumentElement($flowcell_node);

    #generate the XML string
    return $doc->toString(1);

}
