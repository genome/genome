package Genome::Model::Comparison;

use strict;
use warnings;

use Genome;

class Genome::Model::Comparison {
    is => 'Genome::Model',
    has => [
        from_models => {
            is => 'Genome::Model',
            is_input => 1,
            is_many => 1,
            doc => 'the models built in a prior way which are the subject of testing',
        },
        changes => {
            is_many => 1,
            is_input => 1, # is_param => 1,
            is => 'Text', #'UR::BoolExpr' 
            doc => 'changes to the "from" models which are being tested',
        },
    ],
    doc => "pipeline to compare models across processing changes (IN DEVELOPMENT)",
};

sub define_by { return 'Genome::Model::Command::Define::BaseMinimal'; }

sub _help_synopsis {
    return <<EOS;
 genome model define comparison \
    --from id:2890260793/2890224790 \
    --changes "exome_model=''" \
    --processing-profile name="compare nothing" \
    --name my-comparison1 
EOS
}

sub create {
    my $class = shift;
    my $subject = Genome::Taxon->get(name => "unknown");
    my $bx = $class->define_boolexpr(@_)->add_filter(subject => $subject);
    return $class->SUPER::create($bx);
}

sub _execute_build {
    my $self = shift;
    my $build = shift;

    my $from_group_name = $self->name . '.from';
    my @from_models = sort $build->from_models;

    my $from_group = Genome::ModelGroup->get(name => $from_group_name);
    if ($from_group) {
        my @members = sort $from_group->members;
        unless ("@members" eq "@from_models") {
            $self->_rename_model_group($from_group);
            $from_group = undef;
        }
    }
    unless ($from_group) {
        $from_group = Genome::ModelGroup->create(name => $from_group_name);
        for my $from_model (@from_models) {
            $from_group->add_member($from_model);
        }
    }

    my @changes = $build->changes;

    # TODO:
    # pull the logic from this copy command and call it on each $from_model
    # then make an entity representing a pair (UR::Value::Pair?)
    # and set these on the build in some way.
    # Right now we rely on a call to ->members to always sort the same.
    my $to_group_name = $self->name . '.to.' . $self->id;
    my $copy_result = Genome::ModelGroup::Command::Copy->execute(
        from => $from_group,
        to => $to_group_name,
        changes => \@changes,
    );
    my $to_group = Genome::ModelGroup->get(name => $to_group_name);
    unless ($to_group) {
        die $self->error_message("Failed to create model group $to_group_name!");
    }
    my @to_models = $to_group->members;

    for (my $n = 0; $n < $#to_models; $n++) {
        my $from_model = $from_models[$n];
        my $to_model = $to_models[$n];
        my $from_build = $from_model->last_complete_build;
        my $to_build = $to_model->last_complete_build;
        $self->status_message("Compare build " . $from_build->__display_name__ . " to " . $to_build->__display_name__);
    }

    return 1;
}

sub _rename_model_group {
    my ($self, $group) = @_;
    Carp::confess("TODO: implement me!");
}

1;

