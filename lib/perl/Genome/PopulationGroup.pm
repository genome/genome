package Genome::PopulationGroup;

use strict;
use warnings;

use Genome;

class Genome::PopulationGroup {
    is => ['Genome::SampleSource','Genome::Searchable'],
    has => [
        subject_type => { 
            is_constant => 1, 
            is_classwide => 1,
            value => 'sample_group',
        },
        member_hash => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'member_hash' ],
            is_mutable => 1,
            doc => 'MD5 hash created from member IDs, useful for finding population groups by membership',
        },
    ],
    has_many_optional => [
#       This was once on the Subject base class but this caused UR to not know how to join over to
#       SubjectAttributes.
        common_name => {
            calculate_from => 'name',
            calculate => q{ return $name },
        },
        member_ids => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'member' ],
            is_mutable => 1,
            doc => 'IDs of individuals that compose this group',
        },
        members => { 
            is => 'Genome::Individual',
            via => 'attributes',
            to => '_individual',
            where => [ attribute_label => 'member' ],
            doc => 'Genome::Individual objects that compose this group',
        },
        sample_names => {
            via => 'samples',
            to => 'name',
            doc => 'Names of samples that come from individuals in this group',
        },
    ],
    has_optional => [
        taxon_id => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'taxon_id' ],
            is_mutable => 1,
            doc => 'Optional taxon id for this group',
        },
        taxon => { 
            is => 'Genome::Taxon', 
            id_by => 'taxon_id', 
            doc => 'Genome::Taxon object that is assumed to be the taxon for all individuals in this group',
        },
        species_name => { via => 'taxon' },
    ],
    doc => 'A possibly arbitrary grouping of individuals',
};

sub get_source {
    my $self = shift;
    return $self->taxon;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    my @members = $self->members;
    my $member_hash = $self->generate_hash_for_individuals(@members);
    $self->member_hash($member_hash);
    return $self;
}

# Returns all samples derived from the individuals of this group.
# Defining this as a UR property in the class definition above didn't work correctly.
sub samples {
    my $self = shift;
    my @members = $self->members();
    my @samples = Genome::Sample->get(source_id => [ map { $_->id } @members ]);
    return @samples;
}

# Return any other population groups that have the same member hash as the one provided
sub existing_population_groups_with_hash {
    my ($self, $hash) = @_;
    my %params = (member_hash => $hash);
    $params{'id ne'} = $self->id if ref $self; # Exclude self if called as object method
    return Genome::PopulationGroup->get(%params);
}

# Generate an md5 hash based on the ids of the provided individuals, can be called as either an object or class method
sub generate_hash_for_individuals {
    my ($self, @individuals) = @_;
    @individuals = $self->remove_non_unique_individuals(@individuals);
    my @ids = sort map { $_->id } @individuals;
    my $hash = Digest::MD5::md5_hex(@ids); # Returns a valid answer even if individuals is undef
    return $hash;
}

# Filter the provided list of individuals to include only unique individuals
sub remove_non_unique_individuals {
    my ($self, @individuals) = @_;
    return unless @individuals;
    my %unique;
    map { $unique{$_->id} = $_ } @individuals;
    return values %unique;
}

# Returns those members that don't have a matching individual in the provided list
sub find_unmatched_members {
    my ($self, @individuals) = @_;
    return $self->members unless @individuals;

    my %members;
    map { $members{$_->id} = $_ } $self->members;

    for my $individual (@individuals) {
        delete $members{$individual->id} if exists $members{$individual->id};
    }
    return values %members;
}

# Remove individuals from the list that already belong to this group
sub remove_existing_members {
    my ($self, @individuals) = @_;
    return unless @individuals;
    my %individuals;
    map { $individuals{$_->id} = $_ } @individuals;
    for my $member_id ($self->member_ids) {
        delete $individuals{$member_id} if exists $individuals{$member_id};
    }
    return values %individuals;
}

sub add_member {
    my ($self, $individual) = @_;
    return $self->add_members($individual);
}

# Filters the given list of individuals and adds any that pass filtering to the group
sub add_members {
    my ($self, @individuals) = @_;
    my @addable_individuals = $self->remove_non_unique_individuals(@individuals);
    @addable_individuals = $self->remove_existing_members(@addable_individuals);
    return 1 unless @addable_individuals; # If all of the provided individuals are already added/redundant, just do nothing

    for my $addable (@addable_individuals) {
        my $attribute = Genome::SubjectAttribute->create(
            attribute_label => 'member',
            attribute_value => $addable->id,
            subject_id => $self->id,
        );
        unless ($attribute) {
            Carp::confess "Could not add individual " . $addable->id . " to population group " . $self->id;
        }
    }

    my $member_hash = $self->generate_hash_for_individuals($self->members);
    $self->member_hash($member_hash);
    return 1;
}

sub remove_member {
    my ($self, $individual) = @_;
    return $self->remove_members($individual);
}

# Removes members from the group
sub remove_members {
    my ($self, @individuals) = @_;
    my @removed;
    for my $individual (@individuals) {
        my $attribute = Genome::SubjectAttribute->get(
            attribute_label => 'member',
            attribute_value => $individual->id,
            subject_id => $self->id,
        );
        next unless $attribute;
        $attribute->delete;
        push @removed, $individual->id;
    }
    return 1 unless @removed; # Not removing anything is not an error, just return

    my $member_hash = $self->generate_hash_for_individuals($self->members);
    $self->member_hash($member_hash);
    return 1;
}

# Make the groups membership match the provided list of individuals. Assumes that adding/removing
# members from the group is cheap. If this is false, will need to smartly add/remove as needed.
sub change_group_membership {
    my ($self, @individuals) = @_;
    return unless @individuals;
    @individuals = $self->remove_non_unique_individuals(@individuals);
    
    unless ($self->remove_members($self->members)) {
        Carp::confess "Failed to remove members from population group " . $self->__display_name__;
    }

    unless ($self->add_members(@individuals)) {
        Carp::confess "Failed to add individuals " . join(' ', map { $_->__display_name__ } @individuals) .
            " to population group " . $self->__display_name__;
    }

    return 1;
}

1;

