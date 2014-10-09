package Genome::Model::MetagenomicComposition16s::Command; 

use strict;
use warnings;

use Genome;

use Regexp::Common;

class Genome::Model::MetagenomicComposition16s::Command {
    is => 'Genome::Command::Base',
    has => [
        models => {
            is => 'Genome::Model::MetagenomicComposition16s',
            is_many => 1,
            is_optional => 1,
            require_user_verify => 0,
            doc => 'Use the last complete build from these models, resolved via text string',
        },
        builds => {
            is => 'Genome::Model::Build::MetagenomicComposition16s',
            is_many => 1,
            is_optional => 1,
            require_user_verify => 0,
            doc => 'Use these builds, resolved via text string',
        },
    ],
    doc => 'operate on metagenomic-composition-16s models',
};

sub sub_command_category { 'type specific' }

sub _command_name_brief {
    return 'metagenomic-composition-16s';
}

sub _builds {
    my $self = shift;

    my %builds;
    my @models = $self->models;
    for my $model ( @models ) {
        my $build = $model->last_complete_build;
        if ( not $build ) {
            $self->status_message('No last complete build for '.$model->__display_name__);
            next;
        }
        $builds{ $build->id } = $build;
    }

    my @builds = $self->builds;
    for my $build ( @builds ) {
        $builds{ $build->id } = $build;
    }

    if ( not %builds ) {
        $self->error_message('No builds found');
        return;
    }

    return values %builds;
}

1;
