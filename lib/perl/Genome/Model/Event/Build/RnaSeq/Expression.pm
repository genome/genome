package Genome::Model::Event::Build::RnaSeq::Expression;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::RnaSeq::Expression {
    is => ['Genome::Model::Event'],
    is_abstract => 1,
};

sub command_subclassing_model_property {
    return 'expression_name';
}

sub should_bsub { die 'this method shouldn\'t be called.  fix and remove' }
  
1;

