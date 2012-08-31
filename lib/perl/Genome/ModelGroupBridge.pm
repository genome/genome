package Genome::ModelGroupBridge;

use strict;
use warnings;

use Genome;

class Genome::ModelGroupBridge {
    type_name  => 'genome model group bridge',
    table_name => 'GENOME_MODEL_GROUP',
    er_role    => 'bridge',
    id_by      => [
        model_group_id => { is => 'NUMBER', len => 11 },
        model_id       => { is => 'NUMBER', len => 11 },
    ],
    has => [
        model => {
            is              => 'Genome::Model',
            id_by           => 'model_id',
            constraint_name => 'GMG_GM_FK'
        },
        model_group => {
            is              => 'Genome::ModelGroup',
            id_by           => 'model_group_id',
            constraint_name => 'GMG_MG_FK'
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
