package Genome::Test::Factory::SoftwareResult::User;

use strict;
use warnings;

use Genome;
use Genome::Test::Factory::Model::SomaticValidation;
use Genome::Test::Factory::Build;

sub setup_user_hash {
    my $self = shift;
    my @model_params = @_;

    my $model = Genome::Test::Factory::Model::SomaticValidation->setup_object(@model_params);
    my $build = Genome::Test::Factory::Build->setup_object(model_id => $model->id);

    my $user = Genome::Sys->current_user;

    return {
        requestor => $build,
        sponsor   => $user,
    };
}

1;
