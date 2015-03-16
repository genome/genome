package Genome::WorkflowBuilder::StreamGraph;

use strict;
use warnings;
use Genome;
use IPC::Run qw(run);

class Genome::WorkflowBuilder::StreamGraph {
    has => [
        processes => {
            is => 'ARRAY',
            default => [],
        },
        inner_links => {
            is => 'ARRAY',
            default => [],
        },
    ],
};

sub add_process {
    my $self = shift;
    my $process = shift;

    push @{$self->processes}, $process;
}

sub add_link {
    my $self = shift;
    my $link = shift;
    push @{$self->inner_links}, $link;
}

sub get_xml {
    my $self = shift;
    my $doc = XML::LibXML::Document->new();
    $doc->setDocumentElement($self->get_xml_element);

    return $doc->toString(1);
}

sub get_xml_element {
    my $self = shift;
    my $element = XML::LibXML::Element->new('streamgraph');
    for my $process (@{$self->processes}) {
        $element->addChild($process->get_xml_element);
        for my $link ($process->get_link_xml_elements) {
            $element->addChild($link);
        }
    }
    for my $link (@{$self->inner_links}) {
        $element->addChild($link->get_xml_element);
    }
    return $element;
}

sub execute {
    my $self = shift;
    my $xml_file = Genome::Sys->create_temp_file_path;
    Genome::Sys->write_file($xml_file, $self->get_xml);
    run(qw(/home/archive/streamgraph/build/bin/streamgraph run -x), $xml_file, '-o', 'out.xml');#FIXME!
    return 1;
}

1;

