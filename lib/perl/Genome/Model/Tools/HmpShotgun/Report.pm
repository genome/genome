package Genome::Model::Tools::HmpShotgun::Report;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

use Data::Dumper;
use XML::LibXML;

class Genome::Model::Tools::HmpShotgun::Report {
    is  => ['Command'],
    has => [
        working_directory => {
            is => 'String',
            is_input => '1',
            doc => 'The working directory where results will be deposited.',
        },
        delete_intermediates => {
            is => 'Integer',
            is_input =>1,
            is_optional =>1,
            default=>0,
        },
        align_final_file => {
            is => 'String',
            is_input => 1,
            doc => 'The working directory where results will be deposited.',
        },
        final_file => {
            is => 'String',
            is_output => 1,
            is_optional =>1,
            doc => 'The working directory where results will be deposited.',
        },

    ],
};


sub help_brief {
    'Generate reports for HMP Metagenomic Pipeline';
}

sub help_detail {
    return <<EOS
    Generate reports.
EOS
}


sub execute {
    my $self = shift;
    $self->dump_status_messages(1);
    $self->dump_error_messages(1);
    $self->dump_warning_messages(1);

    my $now = UR::Time->now;
    $self->status_message(">>>Starting Report execute() at $now"); 
    
    $self->final_file($self->align_final_file);
    #my $doc = XML::LibXML->createDocument();
    
    #my $root_node = $doc->createElement("hmp-metagenomic-shotgun-report");
    #my $time = UR::Time->now(); 
    #$root_node->addChild( $doc->createAttribute("generated-at",$time) );
    #$doc->setDocumentElement($root_node);
    
    #my $doc_string = $doc->toString(1);
    
    #resolve this
    #my $report_file = $self->working_directory."\reports\final.xml";
    #my $report = Genome::Sys->open_file_for_writing($report_file); 
    
    #print $report $doc_string;
    #$report->close;
    
    
    $self->status_message("<<<Ending Report execute() at ".UR::Time->now); 
    return 1;
}
1;
