package Genome::Command::Delete;

use strict;
use warnings;

use Genome;
      
use Data::Dumper 'Dumper';

class Genome::Command::Delete {
    is => 'Command::V2',
    is_abstract => 1,
    doc => 'CRUD delete command class.',
};

sub _target_name_pl { Carp::confess('Please use CRUD or implement _target_name_pl in '.$_[0]->class); }
sub _target_name_pl_ub { Carp::confess('Please use CRUD or implement _target_name_pl in '.$_[0]->class); }

sub sub_command_sort_position { .4 };

sub help_brief { return $_[0]->_target_name_pl; }

sub help_detail {
    my $class = shift;
    my $target_name_pl = $class->_target_name_pl;
    return "This command deletes $target_name_pl resolved via text string.";
}

sub execute {
    my $self = shift;

    my $target_name_pl_ub = $self->_target_name_pl_ub;
    my @objects = $self->$target_name_pl_ub;

    $self->status_message('Delete '.$self->_target_name_pl.'...');
    for my $obj ( @objects ) {
        $self->_total_command_count($self->_total_command_count + 1);
        my $transaction = UR::Context::Transaction->begin();
        my $display_name = $self->_display_name_for_value($obj);
        my $deleted = eval{ $obj->delete };
        if ($deleted and $transaction->commit) {
            $self->status_message("$display_name: delete");
        }
        else {
            $self->append_error($display_name, "Failed to delete $display_name");
            $transaction->rollback;
        }
    }
    $self->display_command_summary_report();
    $self->status_message('Delete complete. Commiting...');

    return 1; 
}

1;

