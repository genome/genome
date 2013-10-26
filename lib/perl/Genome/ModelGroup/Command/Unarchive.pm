package Genome::ModelGroup::Command::Unarchive;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::Command::Unarchive {
    is => 'Genome::Disk::Command::Allocation::UnarchiveBase',
    has => [
        model_group => {
            is => 'Genome::ModelGroup',
            shell_args_position => 1,
            doc => 'The last completed builds from this model-group will be unarchived',
        },
    ],
};

sub help_brief {
    return "Unarchive all the 'last completed builds' for a model-group";
}

sub help_detail {
    return help_brief();
}

sub _execute {
    my $self = shift;

    my $cmd = Genome::ModelGroup::Command::GetLastCompletedBuilds->execute(
        model_group => $self->model_group,
    );
    my @builds = $cmd->builds;

    my $unarchive_cmd = Genome::Model::Build::Command::Unarchive->create(
        $self->inherited_params,
        builds => \@builds,
    );
    return $unarchive_cmd->execute();
}

# params are the key->value pairs of property->values
sub inherited_params {
    my $self = shift;

    my %params;
    for my $property ($self->inherited_properties) {
        my $name = $property->property_name;
        my $value;
        if ($property->is_many) {
            my @values = $self->$name;
            $value = \@values;
        } else {
            $value = $self->$name;
        }
        $params{$name} = $value;
    }
    return %params;
}

sub inherited_properties {
    my $self = shift;
    my ($parent_class) = @{$self->class->__meta__->{is}}; # is there a better way to get your superclass?
    my @properties = $parent_class->__meta__->properties;

    my @keepers = grep {$_->class_name eq $parent_class} @properties;
    return @keepers;
}

