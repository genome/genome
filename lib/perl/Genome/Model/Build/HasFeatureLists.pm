package Genome::Model::Build::HasFeatureLists;

use strict;
use warnings FATAL => 'all';
use Params::Validate qw(validate validate_pos :types);

class Genome::Model::Build::HasFeatureLists {
    is_abstract => 1,
};

sub feature_list_lookups {
    return {
        target_region => 'get_target_region_feature_list',
        segmental_duplications => 'get_segmental_dupications_feature_list',
        self_chain => 'get_self_chain_feature_list',
    };
}

sub get_feature_list {
    my ($self, $name) = validate_pos(@_, 1, 1);

    my $accessor = $self->feature_list_lookups->{$name};
    unless ($accessor) {
        die $self->error_message("No accessor for name (%s), accessors available for names: %s",
            $name, join(', ', keys %{$self->feature_list_lookups}));
    }

    unless ($self->can($accessor)) {
        die $self->error_message("Accessor (%s) is not defined on this build (%s)",
            $accessor, $self->id);
    }

    my $feature_list = $self->$accessor;
    unless (defined($feature_list) and
        $feature_list->isa('Genome::FeatureList')) {
        die $self->error_message("Couldn't get feature_list named (%s) from accessor (%s)",
            $name, $accessor);
    }

    return $feature_list;
}

sub get_segmental_duplications_feature_list {
    my $self = shift;
    return $self->get_feature_list_from_reference("segmental_duplications");
}

sub get_self_chain_feature_list {
    my $self = shift;
    return $self->get_feature_list_from_reference("self_chain");
}

sub get_feature_list_from_reference {
    die "abstract";
}


