package Genome::Disk::Command::Allocation::Reallocate;

use strict;
use warnings;

use Genome;
use Carp 'confess';

class Genome::Disk::Command::Allocation::Reallocate {
    is => 'Command::V2',
    has => [
        allocations => {
            is => 'Genome::Disk::Allocation',
            doc => 'Allocations to reallocate',
            is_many => 1,
            shell_args_position => 1,
        },
    ],
    has_optional => [
        kilobytes_requested => {
            is => 'Number',
            doc => 'Number of kilobytes that target allocation should reserve, if not ' .
                'provided then the current size of the allocation is used',
        },
        allow_reallocate_with_move => {
            is => 'Boolean',
            default => 0,
            doc => 'Allow the allocation to be moved to a new volume if current volume is too small.',
        },
        grow_only => {
            is => 'Boolean',
            default => 0,
            doc => 'Do not reduce allocation size, only increase it if needed (by du).',
        },
    ],
    doc => 'changes the requested kilobytes for a target allocation',
};

sub help_brief {
    return 'updates kilobytes requested for given allocations';
}

sub help_synopsis { 
    return 'changes the requested kilobytes field on the target allocation';
}

sub help_detail {
    return <<EOS
Changes the requested kilobytes field on the target allocation. If no value
is supplied to this command, the field is set to the current size of the
allocation.
EOS
}

sub execute {
    my $self = shift;
    my @allocations = $self->allocations;
    
    my @errors;
    for my $allocation (@allocations) {
        my $dname = $allocation->__display_name__;
        my %params;
        if ($allocation->is_archived) {
            $self->status_message("Skipping archived allocation $dname");
            next;
        }
        $params{id} = $allocation->id;
        $params{kilobytes_requested} = $self->kilobytes_requested if defined $self->kilobytes_requested;
        $params{allow_reallocate_with_move} = $self->allow_reallocate_with_move;
        $params{grow_only} = $self->grow_only;

        my $transaction = UR::Context::Transaction->begin();
        my $successful = eval {Genome::Disk::Allocation->reallocate(%params) };
        
        if ($successful and $transaction->commit) {
            $self->status_message("Successfully reallocated ($dname).");
        }
        else {
            $self->error_message("Failed to reallocate $dname: $@");
            push @errors, $allocation->__display_name__;
            $transaction->rollback();
        }
    }

    if (@errors) {
        $self->error_message("Failed to reallocate the following allocations: " . join(', ', @errors));
        return 0;
    }
    return 1;
}

1;
