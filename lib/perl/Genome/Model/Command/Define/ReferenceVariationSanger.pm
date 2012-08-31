package Genome::Model::Command::Define::ReferenceVariationSanger;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Define::ReferenceVariationSanger {
    is => 'Genome::Model::Command::Define::HelperDeprecated',
    has => [
        ace_fof =>
        {
	    type => 'String',
            doc => "The full path and filename of the ace_fof for analysis."
        }
    ],

};

sub help_synopsis {
    return <<"EOS"
genome model define reference-variation-sanger -h
EOS
}

sub help_detail {
    return <<"EOS"
This defines the inputs for a new reference-variation-sanger genome model
EOS
}


sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    return $self;
}

sub execute {
    my $self = shift;
    unless(defined $self->ace_fof) {
	$self->error_message("ace_fof is required to define model");
	return;
    }
    my $ace_fof = $self->ace_fof;
    unless (-f $ace_fof) {$self->error_message("couldn't find your ace_fof");return;}

    my $super = $self->super_can('_execute_body');
    $super->($self,@_);    

    my $model = Genome::Model->get($self->result_model_id);
    $model->add_input(name => 'ace_fof', value_class_name => 'UR::Value', value_id => $ace_fof);

    return 1;
 
}

1;
