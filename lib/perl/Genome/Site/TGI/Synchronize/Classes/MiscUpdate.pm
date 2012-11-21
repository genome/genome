package Genome::Site::TGI::Synchronize::Classes::MiscUpdate; 

use strict;
use warnings;

=pod
SUBJECT_CLASS_NAME    VARCHAR2 (255)                   {null} NOT NULL pk
SUBJECT_ID            VARCHAR2 (255)                   {null} NOT NULL pk
SUBJECT_PROPERTY_NAME VARCHAR2 (255)                   {null} NOT NULL pk
EDIT_DATE             TIMESTAMP(6)                     {null} NOT NULL pk
OLD_VALUE             VARCHAR2 (1000)                  {null} {null}   
NEW_VALUE             VARCHAR2 (1000)                  {null} {null}   
EDITOR_ID             VARCHAR2 (255)                   {null} NOT NULL 
DESCRIPTION           VARCHAR2 (255)                   {null} NOT NULL 
IS_RECONCILED         NUMBER   (1)                     0      NOT NULL 
=cut

class Genome::Site::TGI::Synchronize::Classes::MiscUpdate {
    is => 'UR::Object',
    table_name => 'GSC.misc_update',
    id_by => [
        subject_class_name => {
            is => 'Text', 
            doc => 'the table in which the change occurred' 
        },
        subject_id => { 
            is => 'Text', 
            doc => 'the primary key of the row that changed' 
        },
        subject_property_name => {
            is => 'Text', 
            doc => 'the column whose value changed' 
        },
        edit_date => { 
            is => 'Date', 
            doc => 'the time of the change' 
        },
    ],
    has => [
        editor_id => { 
            is => 'Text', 
            doc => 'the unix account that made the change' 
        },
        description => { 
            is => 'Text', 
            valid_values => [qw/ INSERT UPDATE DELETE /], 
            doc => 'the type of change' 
        },
        is_reconciled => { 
            is => 'Boolean', 
            doc => 'Indicates if the update has been applied to our tables'
        },
    ],
    has_optional => [
        old_value => { 
            is => 'Text', 
            doc => 'the value which was changed' 
        },
        new_value => { 
            is => 'Text', 
            doc => 'the value to which old_value was changed' 
        },
    ],
    data_source => 'Genome::DataSource::GMSchema',
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

