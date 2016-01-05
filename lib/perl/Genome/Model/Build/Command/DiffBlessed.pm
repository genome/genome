package Genome::Model::Build::Command::DiffBlessed;

use strict;
use warnings;

use Genome;
use File::Spec;
use YAML;

class Genome::Model::Build::Command::DiffBlessed {
    is => 'Genome::Interfaces::Comparable::Command::DiffBlessed',
    has => [
        new_build => {
            is => 'Genome::Model::Build',
        },
        perl_version => {
            is => 'Text',
            valid_values => ['5.8', '5.10'],
            default_value => '5.10',
        },
    ],
};

sub blessed_object {
    my $self = shift;
    return $self->blessed_build;
}

sub new_object {
    my $self = shift;
    return $self->new_build;
}

sub blessed_build {
    my $self = shift;
    my $model_name = $self->new_build->model_name;
    return $self->retrieve_blessed_build($model_name);
}

sub retrieve_blessed_build {
    my ($self, $model_name) = @_;
    my $blessed_id = $self->blessed_id($model_name);
    my $blessed_build = Genome::Model::Build->get(
        id => $blessed_id,);
    return $blessed_build;
}

