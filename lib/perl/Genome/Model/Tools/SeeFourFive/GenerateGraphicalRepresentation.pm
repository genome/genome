package Genome::Model::Tools::SeeFourFive::GenerateGraphicalRepresentation;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Workflow;
use FileHandle;

class Genome::Model::Tools::SeeFourFive::GenerateGraphicalRepresentation{
    is => 'Command',
    has => [
        c_file =>
        {
            is => 'String',
            doc => 'c4.5 tree output'
        },
        ]
};

sub help_synopsis {
    "creates a png graph of the passed decision tree"
}

sub execute {
    my $self=shift;
    
    my $c45_object = Genome::Model::Tools::SeeFourFive::Tree ->create();
    $c45_object->c45_file($self->c_file);
    $c45_object->load_trees;
    my $graph = $c45_object->as_graphviz_obj;
    
    #print the graph to temp
    if($graph) {
        $graph->as_png("/tmp/c4.5graph.png");
    }
    else {
        return;
    }
    
    return 1;
}

1;
