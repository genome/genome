package Genome::Disk::Command::Allocation::Deallocate;

use strict;
use warnings;

use Genome;
use Carp 'confess';

class Genome::Disk::Command::Allocation::Deallocate {
    is => 'Command::V2',
    has_many => [        
        allocations => {
            is => 'Genome::Disk::Allocation',
            shell_args_position => 1,
            doc => 'Allocations to delete',
        },
    ],
    doc => 'removes target allocation and deletes its directories',
};

sub help_brief {
    return 'removes the target allocation and deletes its directories';
}

sub help_synopsis {
    return 'removes the target allocation and deletes its directories';
}

sub help_detail {
    return 'removes the target allocation and deletes its directories';
}

sub execute { 
    my $self = shift;
    my @allocations = $self->allocations;

    my @errors;
    for my $allocation (@allocations) {
        my $display_name = $allocation->__display_name__;
        my $transaction = UR::Context::Transaction->begin();

        my $successful = Genome::Disk::Allocation->delete(id => $allocation->id);

        if ($successful and $transaction->commit) {
            $self->status_message("Successfully deallocated ($display_name)");
        }
        else {
            $self->error_message("Failed to deallocate ($display_name): $@");
            push @errors, $display_name;
            $transaction->rollback;
        }
    }

    if (@errors) {
        $self->error_message("Failed to deallocate the following allocations: " . join(', ', @errors));
        return 0;
    }
    return 1;
}
    
1;
