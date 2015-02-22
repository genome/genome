package Genome::Model::ReferenceAlignment::Command::ExampleOneModel; 

class Genome::Model::ReferenceAlignment::Command::ExampleOneModel {
    is => 'Command::V2',
    has => [
        model => { shell_args_position => 1, is => 'Genome::Model', id_by => 'model_id' }, 
    ],
    doc => 'example command which works on one model'
};

sub execute {
    my $self = shift;
    my $model = $self->model;
    print "model $model\n";
    return 1;
}

1;
