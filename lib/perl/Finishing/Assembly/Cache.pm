package Finishing::Assembly::Cache;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;

my %data    :name(data:o) :ds(hashref) :default({}) :empty_ok(1);

sub ids
{
    my $self = shift;

    return keys %{ $self->data };
}

sub id_exists
{
    my ($self, $id) = @_;

    $self->fatal_msg
    (
        "No id to test exists", { caller_level => 1 }
    ) unless $id;
    
    return exists $self->data->{$id};
}

sub get
{
    my ($self, $id, $attr) = @_;

    $self->fatal_msg
    (
        "No id to get attribute", { caller_level => 1 }
    ) unless defined $id;
    
    $self->fatal_msg
    (
        "No attribute to get for id ($id)", { caller_level => 1 }
    ) unless defined $attr;
 
    return $self->data->{$id}->{$attr};
}

sub get_all_properties_for_id
{
    my ($self, $id) = @_;

    $self->fatal_msg
    (
        "No id to get attribute", { caller_level => 1 }
    ) unless defined $id;
    
    # TODO error?
    return unless exists $self->data->{$id};

    return %{ $self->data->{$id} };
}

sub set
{
    my ($self, $id, $attr, $val) = @_;
    
    $self->fatal_msg("No id to set attribute", { caller_level => 1 }) unless defined $id;
    $self->fatal_msg("No attribute to set for id ($id)", { caller_level => 1 }) unless defined $attr;
    $self->fatal_msg
    (
        "No value to set attribute ($attr) for id ($id)", { caller_level => 1 }
    ) unless defined $val;
    
    $self->data->{$id}->{$attr} = $val;
    
    return $val;
}

sub set_all_properties_for_id
{
    my ($self, $id, $props) = @_;

    $self->fatal_msg
    (
        "No id to get attribute", { caller_level => 1 }
    ) unless defined $id;
    $self->fatal_msg("No propertires to set for id ($id)", { caller_level => 1 }) unless $props;
    $self->fatal_msg("Propertires to set needs to be a hash ref for id ($id)", { caller_level => 1 }) unless %$props;

    $self->data->{$id} = \%{ %$props };

    return 1;
}

sub push
{
    my ($self, $id, $attr, @vals) = @_;
    
    $self->fatal_msg("No id to push data", { caller_level => 1 }) unless defined $id;
    $self->fatal_msg("No attribute to push for id ($id)", { caller_level => 1 }) unless defined $attr;
    $self->fatal_msg
    (
        "No values to push to attribute ($attr) for id ($id)", { caller_level => 1 }
    ) unless @vals;

    if ( my $stored_val = $self->data->{$id}->{$attr} )
    {
        my $ref = ref($stored_val);
        if ( $ref and $ref ne 'ARRAY' )
        {
            $self->data->{$id}->{$attr} = [ $stored_val ];
        }
    }

    push @{ $self->data->{$id}->{$attr} } , @vals;

    return $self->data->{$id}->{$attr};
}

sub move
{
    my ($self, $id, $new_id) = @_;

    $self->fatal_msg
    (
        "No id to move data", { caller_level => 1 }
    ) unless defined $id;

    $self->fatal_msg
    (
        "Can't find data for id ($id)", { caller_level => 2 }
    ) unless exists $self->data->{$id};

    $self->fatal_msg
    (
        "No new id to move data for id ($id)", { caller_level => 2 }
    ) unless defined $new_id;

    $self->fatal_msg
    (
        "New id ($new_id) already exists, can't move data to it", { caller_level => 2 }
    ) if exists $self->data->{$new_id};
    
    return $self->data->{$new_id} = delete $self->data->{$id};
}

sub delete
{
    my ($self, $id) = @_;

    $self->fatal_msg
    (
        "No id to delete data", { caller_level => 1 }
    ) unless defined $id;

    $self->fatal_msg # needed?
    (
        "Can't find data for id ($id) to delete", { caller_level => 1 }
    ) unless exists $self->data->{$id};
    print "deleting id: $id attrs: ",keys %{$self->data->{$id}},"\n";
    return delete $self->data->{$id};
}

sub undef_attr
{
    my ($self, $id, $attr) = @_;
    
    $self->fatal_msg("No id to set attribute", { caller_level => 1 }) unless defined $id;
    $self->fatal_msg("No attribute to set for id ($id)", { caller_level => 1 }) unless defined $attr;
    
    $self->data->{$id}->{$attr} = undef;
    
    return 1;
}

1;

#$HeadURL$
#$Id$
