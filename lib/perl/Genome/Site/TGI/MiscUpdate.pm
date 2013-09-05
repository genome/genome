package Genome::Site::TGI::MiscUpdate; 

use strict;
use warnings;

class Genome::Site::TGI::MiscUpdate {
    table_name => 'misc_update',
    id_by => [
        table_name => {
            is => 'TExt',
            column_name => 'SUBJECT_CLASS_NAME',
            is => 'Text', 
            len => 255, 
            doc => 'the table in which the change occurred' 
        },
        updated_row_primary_key => { 
            is => 'Text', 
            column_name => 'SUBJECT_ID',
            len => 255, 
            doc => 'the primary key of the row that changed' 
        },
        edit_date => { 
            is => 'Date', 
            doc => 'the time of the change' 
        },
    ],
    has => [
        editor_id => { 
            is => 'Text', 
            len => 255, 
            doc => 'the unix account that made the change' 
        },
        column_name => {
            is => 'Text', 
            column_name => 'SUBJECT_PROPERTY_NAME',
            len => 255, 
            doc => 'the column whose value changed' 
        },
        type => { 
            is => 'Text', 
            column_name => 'DESCRIPTION',
            len => 255, 
            valid_values => ['INSERT', 'UPDATE', 'DELETE'], 
            doc => 'the type of change (we do not currently track inserts)' 
        },
        is_reconciled => { 
            is => 'Boolean', 
            default => 0, 
            doc => 'Indicates if the update has been applied to our tables'
        },
    ],
    has_optional => [
        old_value => { 
            is => 'Text', 
            len => 1000, 
            doc => 'the value which was changed' 
        },
        new_value => { 
            is => 'Text', 
            len => 1000, 
            doc => 'the value to which old_value was changed' 
        },
    ],
    data_source => 'Genome::DataSource::Dwrac',
    doc => 'The MISC_UPDATE table tracks changes to certain other tables in the gsc schema.'
};

sub __display_name__ {
    my $self = shift;
    my $table = $self->table_name;
    my $column = $self->column_name || 'unknown_column';
    my $old_value = $self->old_value || 'null';
    my $new_value = $self->new_value || 'null';
    return "$table.$column($old_value => $new_value)";
}

1;

