package Genome::Model::Tools::Tcga::AppendColumns;

use strict;
use warnings;

use Genome;
use Command;
use FileHandle;

class Genome::Model::Tools::Tcga::AppendColumns {
    is => 'Command',
    has => [
    f_of_columns => {
        type => 'String',
        is_optional => 0,
        doc => "File of chromosome and position and of columns of the query",
    },
    f_to_add_columns => {
        type => 'String',
        is_optional => 0,
        doc => 'The files to query and be appended',
    },
    columns => {
        type => 'String',
        is_optional => 0,
        doc => 'Column numbers of f_of_columns, with comma separated',
    },
    ],
    has_optional => [
    positions_a => {
        type => 'String',
        is_optional => 1,
        default => '0,1,2',
        doc => 'Column numbers of chromosome, start and end of file of the columns',
    },        
    positions_b => {
        type => 'String',
        is_optional => 1,
        default => '0,1,2',
        doc => 'Column numbers of chromosome, start and end of file to add the columns',
    },
    ],
};

sub help_brief {
}

sub help_detail {
}

sub execute {
    my $self=shift;
    my %query_positions;

    my @pos = split ",", $self->positions_a;
    my @col = split ",", $self->columns;

    #read in the query positions
    my $handle = new FileHandle;
    $handle->open($self->f_of_columns,"r");
    unless($handle) {
        $self->error_message("Unable to open file with the columns: " . $self->f_of_columns);
        return;
    }

    while(my $line = $handle->getline) {
        chomp $line;
        my @all_columns_a = split "\t", $line;
        my $chromosome = $all_columns_a[$pos[0]];
        my $start = $all_columns_a[$pos[1]];
        my $stop = $all_columns_a[$pos[2]];
        my $extra_column = "";
        for(my $i = 0; $i <= $#col; $i ++){
            $extra_column = $extra_column . $all_columns_a[$col[$i]];
            $extra_column = $extra_column . ",," if($i != $#col);
        }
        $query_positions{$chromosome}{$start}{$stop} = $extra_column;
    }
    $handle->close;

    my $handle_b = new FileHandle;
    $handle_b->open($self->f_to_add_columns,"r");
    unless($handle_b){
        $self->error_message("Unable to open file with the columns to be added" . $self->f_to_add_columns);
        return;
    }

    my @pos_b = split ",", $self->positions_b;
    while(my $line = $handle_b->getline) {
        chomp $line;
        my @all_columns_a = split "\t", $line;
        my $chromosome = $all_columns_a[$pos_b[0]];
        my $start = $all_columns_a[$pos_b[1]];
        my $stop = $all_columns_a[$pos_b[2]];
        print $line;
        if(defined $query_positions{$chromosome}{$start}{$stop}){
            my $extra_column_to_add = $query_positions{$chromosome}{$start}{$stop};
            my @extra_column_to_add_ = split ",,", $query_positions{$chromosome}{$start}{$stop};
            for(my $i = 0; $i <= $#extra_column_to_add_; $i++){
                print "\t". $extra_column_to_add_[$i];
            }
        }
        print "\n";
    }

    $handle_b->close;

    return 1;
}


1;




