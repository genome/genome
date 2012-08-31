package Genome::Model::Build::ExecuteBuildWrapper;

use strict;
use warnings;
use Genome;

# this command is not intended to be executed by users
# the only purpose of this command is to call _execute_build, or some other method, on the model
# so it can be turned into a one-step workflow.

class Genome::Model::Build::ExecuteBuildWrapper {
    is  => 'Command',
    has_input => [
        build_id => {
            is  => 'Number',
            doc => 'specify the build by id'
        },
        method_name => {
            is => 'Text',
            is_optional => 1,
            default_value => '_execute_build',
            doc => 'the method on the processing profile to call, passing a build (defaults to _execute_build())',
        }
    ]
};

sub bsub_rusage {
    return "-R 'select[model!=Opteron250 && type==LINUX64] rusage[tmp=90000:mem=16000]' -M 16000000";
}

sub execute {
    my $self = shift;


    my $build = Genome::Model::Build->get($self->build_id) 
      or die 'cannot load build object for ' . $self->build_id;

    my $method_name = $self->method_name;

    my $model = $build->model;
    my $pp = $build->processing_profile; #TODO remove after '_execute_build's are all moved to model subclasses.

    my $rv;
    if ($model->can($method_name)) {
        $rv = $model->$method_name($build);
    }
    elsif ($pp->can($method_name)) {
        $rv = $pp->$method_name($build);
    }

    die $method_name . ' returned undef' if !defined $rv;
    die $method_name . ' returned false' if !$rv;

    return $rv;
}

1;

