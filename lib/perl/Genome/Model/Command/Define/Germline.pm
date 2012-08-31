package Genome::Model::Command::Define::Germline;

use strict;
use warnings;
use Genome;

class Genome::Model::Command::Define::Germline {
    is => 'Genome::Model::Command::Define::HelperDeprecated', 
    has => [
        source_model_id => {
            is => 'Text',
        },
        source_model => {
            is => 'Genome::Model',
            id_by => 'source_model_id',
            doc => 'Model that germline model will analyze',
        }
    ],
};

sub type_specific_parameters_for_create {
    my $self = shift;
    my @params = ();
    push @params, (
        source_model_id => $self->source_model_id
    );
    return @params;
}

sub execute {
    my $self = shift;
    my $result =  $self->SUPER::_execute_body(@_);
    return unless $result;

    my $model_id = $self->result_model_id;
    my $model = Genome::Model->get($model_id);
    Carp::confess "Could not retrieve model $model_id!" unless $model;

    my $model_subject = $model->subject;
    my $source_model_subject = $model->source_model->subject;

    unless ($model_subject->id eq $source_model_subject->id and $model_subject->class eq $source_model_subject->class) {
        Carp::confess 'Provided subject (ID ' . $model_subject->id . ', class ' . $model_subject->class .
            ") does not match source model's subject (ID " . $source_model_subject->id . ', class ' . $source_model_subject->class . ')';
    }

}
1;

