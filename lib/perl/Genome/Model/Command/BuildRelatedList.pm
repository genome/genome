package Genome::Model::Command::BuildRelatedList;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::BuildRelatedList {
    is => 'UR::Object::Command::List',
    is_abstract => 1,
    has_optional => [
        filter => {
            # override shell_args_position from UR::Object::Command::List
            shell_args_position => undef,
        },
        build_spec => {
            is => 'Text',
            shell_args_position => 1,
            doc => "The id of the build, or the id or name of a model, or a filter",
        },
    ],
};

sub create {
    # TODO: pull this up into Genome::Model::Command so it's fast/easy to run model/build centric things.
    # 2011-10: instead can now just use Command::V2 but that isn't integrated into lister yet
    my $class = shift;

    # TODO: get rid of the ' ' key in the construction params
    # It was a hack to get around not having "shell_args_position" in the properties.
    # 2011-10: Now shell_args_position is set but this does more, e.g. match by name or ID
    # lots of overlap with Command::V2 but that isn't integrated into lister yet
    my %params = @_;
    delete $params{' '};
    my $bx = $class->define_boolexpr(%params);

    if (my $build_spec = $bx->value_for("build_spec")) {
        my $filter_updated;
        my $filter_prefix = ($bx->value_for("filter") ? ($bx->value_for("filter") . ',') : '');
        if ($build_spec !~ /\D/) {
            my $build = Genome::Model::Build->get($build_spec);
            if ($build) {
                $filter_updated = $filter_prefix . 'build_id=' . $build->id
            }
            else {
                my $model = Genome::Model->get($build_spec);
                if ($model) {
                    $filter_updated = $filter_prefix . 'model_id=' . $model->id;
                }
            }
        }
        unless ($filter_updated) {
            my @models = Genome::Model->get(name => $build_spec);
            unless (@models) {
                @models = Genome::Model->get("name like" => $build_spec . '%');
            }
            if (@models == 1) {
                $filter_updated = $filter_prefix . 'model_id=' . $models[0]->id;
            }
            elsif (@models > 1) {
                $filter_updated = $filter_prefix . 'model_name~' . $build_spec . '%';
            }
            else {
                $filter_updated = $filter_prefix . $build_spec;
            }
        }
        if ($filter_updated) {
            return $class->SUPER::create($bx->params_list, filter => $filter_updated);
        }
    }
    return $class->SUPER::create(@_);
}

1;

