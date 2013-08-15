package Genome::ModelGroupBridge;

use strict;
use warnings;

use Genome;

class Genome::ModelGroupBridge {
    table_name => 'model.model_group_bridge',
    er_role => 'bridge',
    type_name => 'genome model group bridge',
    id_by => [
        model_group_id => {
            is => 'Text',
            len => 64,
            is_deprecated => 1,
        },
        model_id => {
            is => 'NUMBER',
            len => 11,
        },
    ],
    has => [
        model => {
            is => 'Genome::Model',
            id_by => 'model_id',
            constraint_name => 'GMG_GM_FK',
        },
        group => {
            is => 'Genome::ModelGroup',
            id_by => 'model_group_id',
            constraint_name => 'GMG_MG_FK',
            is_deprecated => 1,
        },
        model_group => {
            is => 'Genome::ModelGroup',
            id_by => 'model_group_id',
            is_deprecated => 1,
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};


if ($INC{"Genome/Search.pm"}) {
    my $queue_model_group_callback = sub {
        my ($self) = @_;
        my $mg = $self->model_group;
        Genome::Search::Queue->create(
            subject_id => $mg->id,
            subject_class => $mg->class,
        );
    };
    Genome::ModelGroupBridge->add_observer(
        callback => $queue_model_group_callback,
        aspect => 'create',
    );
    Genome::ModelGroupBridge->add_observer(
        callback => $queue_model_group_callback,
        aspect => 'delete',
    );
}

1;
