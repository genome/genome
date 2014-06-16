package Genome::Utility::ObjectWithAllocations;

use strict;
use warnings;

use Genome;

require List::MoreUtils;

class Genome::Utility::ObjectWithAllocations {
    is => 'UR::Object',
    is_abstract => 1,
    has => {
        disk_allocations => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            reverse_as => 'owner',
        },
    },
};

sub delete {
    my $self = shift;
    $self->_create_deallocate_observer;
    return $self->SUPER::delete;
}

sub _create_deallocate_observer {
    my $self = shift;

    my $deallocator;
    $deallocator = sub {
        $self->deallocate;
        UR::Context->cancel_change_subscription(
            'commit', $deallocator
        );
    };
    UR::Context->create_subscription(
        method => 'commit',
        callback => $deallocator
    );

    return 1;
}

sub reallocate {
    # FIXME
    #  what to do 'with_move'?
    #  do we care if reallocate fails?
    my ($self, $with_move) = @_; 

    for my $disk_allocation ( $self->disk_allocations ) {
        my $reallocate_ok = eval { $disk_allocation->reallocate };
        next if $reallocate_ok;
        $self->warning_message($@) if $@;
        $self->warning_message('Continuing, but failed to reallocate disk allocation: '.$disk_allocation->__display_name__);
    }

    return 1;
}

sub deallocate { 
    # FIXME
    #  do we care if deallocate fails?
    my $self = shift;

    for my $disk_allocation ( $self->disk_allocations ) {
        print "$disk_allocation\n";
        next if $disk_allocation->isa('UR::DeletedRef');
        my $deallocate_ok = eval{ $disk_allocation->deallocate; };
        next if $deallocate_ok;
        $self->warning_message($@) if $@;
        $self->warning_message('Continuing, but failed to deallocate disk allocation: '.$disk_allocation->__display_name__);
    }

    return 1;
}

sub archivable {
    my $self = shift;

    my @disk_allocations = $self->disk_allocations;
    return 0 if not @disk_allocations;

    for my $disk_allocation ( @disk_allocations ) {
        return 0 if not $disk_allocation->archivable;
    }

    return 1;
}

sub is_archived {
    my $self = shift;

    my @disk_allocations = $self->disk_allocations;
    return 0 if not @disk_allocations;

    for my $disk_allocation ( @disk_allocations ) {
        return 1 if $disk_allocation->is_archived;
    }

    return 0;
}

sub associated_disk_allocations {
    my $self = shift;

    my @associated_disk_allocations = $self->disk_allocations;
    my @addtional_associated_disk_allocations = $self->_additional_associated_disk_allocations;
    push @associated_disk_allocations, @addtional_associated_disk_allocations if @addtional_associated_disk_allocations;

    return List::MoreUtils::uniq(@associated_disk_allocations);
}
sub _additional_associated_disk_allocations {}

sub unarchive {
    my $self = shift;
    
    my $unarchive_by_command_class = $self->_unarchive_by_command_class;
    return if $unarchive_by_command_class;

    for my $disk_allocation ( $self->disk_allocations ) {
        my $unarchive_ok = eval{ $disk_allocation->unarchive; };
        next if $unarchive_ok;
        $self->error_message($@) if $@;
        $self->error_message('Failed to unarchive disk allocation! '.$disk_allocation->__display_name__);
        return;
    }

    return 1;
}

sub _unarchive_by_command_class {
    my $self = shift;

    # FIXME there has to be a beeter way to get the unarchive command class
    my %classes_and_unarchive_command_classes = (
        'Genome::Model::Build' => 'Genome::Model::Build::Command::Unarchive',
    );
    for my $class ( keys %classes_and_unarchive_command_classes ){
        return $classes_and_unarchive_command_classes{$class} if $self->isa($class);
    }

    return;
}

1;

