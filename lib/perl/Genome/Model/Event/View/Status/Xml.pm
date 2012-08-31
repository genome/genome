package Genome::Model::Event::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has => [
        _doc    => {
            is_transient => 1,
            doc => 'the XML::LibXML document object used to build the content for this view'
        },
    ],
};

sub _generate_content {
    my $self = shift;

    #create the XML doc and add it to the object
    my $doc = XML::LibXML->createDocument();
    $self->_doc($doc);

    my $subject = $self->subject;
    return unless $subject;

    my $event_node = $self->get_event_node($subject);
    $doc->setDocumentElement($event_node);

    #generate the XML string
    return $doc->toString(1);
}

sub get_event_node {
    my $self = shift;
    my $event = shift;
    my $doc = $self->_doc;

    my $event_node = $self->anode("genome-model-event","id",$event->id);
    $event_node->addChild( $doc->createAttribute("command_class",$event->class));
    $event_node->addChild( $self->tnode("event_status",$event->event_status));

    my $lsf_job_id = $event->lsf_job_id;

    my $lsf_job_status = $self->get_lsf_job_status($lsf_job_id);

    $event_node->addChild( $self->tnode("lsf_job_id",$lsf_job_id));
    $event_node->addChild( $self->tnode("lsf_job_status",$lsf_job_status));
    $event_node->addChild( $self->tnode("date_scheduled",$event->date_scheduled));
    $event_node->addChild( $self->tnode("date_completed",$event->date_completed));
    $event_node->addChild( $self->tnode("elapsed_time", $self->calculate_elapsed_time($event->date_scheduled,$event->date_completed) ));
    $event_node->addChild( $self->tnode("instrument_data_id",$event->instrument_data_id));
    my $err_log_file = $event->resolve_log_directory ."/".$event->id.".err";
    my $out_log_file = $event->resolve_log_directory ."/".$event->id.".out";
    $event_node->addChild( $self->tnode("output_log_file",$out_log_file));
    $event_node->addChild( $self->tnode("error_log_file",$err_log_file));

    # get alignment director[y|ies] and filter description
    if($event->instrument_data_id) {
        my $instrument_data = Genome::InstrumentData->get($event->instrument_data_id);
        $event_node->addChild( $self->get_instrument_data_node($instrument_data));

        my @inputs = $event->model->input_for_instrument_data_id($event->instrument_data_id);

        if (scalar @inputs > 0) {
            # find the events with matching instrument_data_ids
            my @adirs;

            my $processing_profile = $event->model->processing_profile;
            for my $input (@inputs) {
                my $instrument_data = $input->value;
                if ($instrument_data->id == $event->instrument_data_id) {
                    my $alignment;
                    my %segment_info;
                    if (defined $event->instrument_data_segment_id) {
                        $segment_info{instrument_data_segment_id} = $event->instrument_data_segment_id; 
                        $segment_info{instrument_data_segment_type} = $event->instrument_data_segment_type; 
                    }
                    eval{ ($alignment) = $processing_profile->results_for_instrument_data_input($input, %segment_info)};

                    if ($@) {
                        chomp($@);
                        push(@adirs, $@);
                    }

                    if (defined($alignment)) {
                        push(@adirs, $alignment->output_dir);

                        # look for a filter description
                        if ($input->filter_desc) {
                            $event_node->addChild( $self->tnode("filter_desc", $input->filter_desc));
                        }
                    }
                }
            }
            # handle multiple alignment directories
            if (scalar @adirs > 1) {
                my $i = 1;
                for my $adir (@adirs) {
                    $event_node->addChild( $self->tnode("alignment_directory_" . $i, $adir));
                    $i++;
                }
            } else {
                $event_node->addChild( $self->tnode("alignment_directory", $adirs[0]));
            }
        }
    }
    return $event_node;
}

sub get_instrument_data_node {

    my $self = shift;
    my $object = shift;

    my $id = $self->anode("instrument_data","id", $object->id);
    for (qw/flow_cell_id project_name run_name run_identifier read_length library_name library_id lane subset_name run_type gerald_directory id/) {
        if ($object->class ne 'Genome::InstrumentData::Imported' && $object->can($_)) {
            $id->addChild($self->tnode($_, $object->$_));
        } else {
            $id->addChild($self->tnode($_, "N/A")); 
        }     
    }

    return $id;

}

#helper methods.  just pass through to the more descriptive names
#anode = attribute node
sub anode {
    my $self = shift;
    return $self->create_node_with_attribute(@_);
}

sub create_node_with_attribute {

    my $self = shift;
    my $node_name = shift;
    my $attr_name = shift;
    my $attr_value = shift;

    my $doc = $self->_doc;

    my $node = $doc->createElement($node_name);
    $node->addChild($doc->createAttribute($attr_name,$attr_value));
    return $node;

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

    my $doc = $self->_doc;

    my $node = $doc->createElement($node_name);
    if ( defined($node_value) ) {
        $node->addChild($doc->createTextNode($node_value));
    }
    return $node;

}

sub add_attribute {
    my $self = shift;
    my $node = shift;
    my $attr_name = shift;
    my $attr_value = shift;

    my $doc = $self->_doc;

    $node->addChild($doc->createAttribute($attr_name,$attr_value) );
    return $node;

}

sub calculate_elapsed_time {
    my $self = shift;
    my $date_scheduled = shift;
    my $date_completed = shift;

    my $diff;

    if ($date_completed) {
        $diff = UR::Time->datetime_to_time($date_completed) - UR::Time->datetime_to_time($date_scheduled);
    } else {
        $diff = time - UR::Time->datetime_to_time( $date_scheduled);
    }

    # convert seconds to days, hours, minutes
    my $seconds = $diff;
    my $days = int($seconds/(24*60*60));
    $seconds -= $days*24*60*60;
    my $hours = int($seconds/(60*60));
    $seconds -= $hours*60*60;
    my $minutes = int($seconds/60);
    $seconds -= $minutes*60;

    my $formatted_time;
    if ($days) {
        $formatted_time = sprintf("%d:%02d:%02d:%02d",$days,$hours,$minutes,$seconds);
    } elsif ($hours) {
        $formatted_time = sprintf("%02d:%02d:%02d",$hours,$minutes,$seconds);
    } elsif ($minutes) {
        $formatted_time = sprintf("%02d:%02d",$minutes,$seconds);
    } else {
        $formatted_time = sprintf("%02d:%02d",$minutes,$seconds);
    }

    return $formatted_time;

}


our $JOB_TO_STATUS;
sub load_lsf_job_status {
    my $self = shift;

    # NOTE: This caches the lsf job data for as long as the process stays alive.
    # For now, this thing is run from a regular CGI script, so the process dies pretty
    # quickly.  But if things change so that the process lives for a while, then
    # a different cache aging mechanism should be set up
    unless ($JOB_TO_STATUS) {
        $JOB_TO_STATUS = {};

        my $lsf_file = '/gsc/var/cache/testsuite/lsf-tmp/bjob_query_result.txt';
        my $fh = IO::File->new($lsf_file);
        my $lsf_file_data = do { local( $/ ) ; <$fh> } ;
        $fh->close;
        while ($lsf_file_data =~ m/^(\S+)\s+(\S+)\s+(\S+).*?\n/gm) {
            $JOB_TO_STATUS->{$1} = $3;
        }
        delete $JOB_TO_STATUS->{'JOBID'};
    }
    return $JOB_TO_STATUS;
}



sub get_lsf_job_status {
    my $self = shift;
    my $lsf_job_id = shift;

    my $result;

    if ( defined($lsf_job_id) ) {

        #check the user specified flag to determine how to retrieve lsf status
#        if ($self->use_lsf_file) {
#            #get the data from the preloaded hash of lsf info (from file)
#            #my %job_to_status = %{$self->_job_to_status};
#            #$result = $job_to_status {$lsf_job_id};
#            $result = $self->{'_job_to_status'}->{$lsf_job_id};
#            if (!defined($result) ) {
#                $result = "UNAVAILABLE";
#            }
#        } else {
            #get the data directly from lsf via bjobs command
            my @lines = `bjobs $lsf_job_id 2>/dev/null`;
            #parse the bjobs output.  get the 3rd field of the 2nd line.
            if ( (scalar(@lines)) > 1) {
                my $line = $lines[1];
                my @fields = split(" ",$line);
                $result = $fields[2];
            } else {
                #if there are no results from bjobs, lsf forgot about the job already.
                $result = "UNAVAILABLE";
            }
#        }

    } else {
        #if the input LSF ID is not defined, mark it as unscheduled.
        $result = "UNSCHEDULED";
    }
    return $result;

    #NOTES:  UNSCHEDULED means that an LSF ID exists, but LSF did not have any status on it.  Probably because it was executed a while ago.
    #        UNAVAILABLE means that an LSF ID does NOT exist.
}

1;
