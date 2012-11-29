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
    has_transient => [
        result => { is => 'Text', },
        status => { is => 'Text', },
    ],
    data_source => 'Genome::DataSource::GMSchema',
};

sub lims_table_name {
    my $self = shift;

    my $subject_class_name = $self->subject_class_name;
    my ($schema, $lims_table_name) = split(/\./, $subject_class_name);
    if ( not $schema ) {
        $self->error_message('Failed to get schema from subject class name => '.$subject_class_name);
        return;
    }
    if ( not $lims_table_name ) {
        $self->error_message('Failed to get lims table name from subject class name => '.$subject_class_name);
        return;
    }

    return $lims_table_name;
}

my %subject_class_names_to_genome_class_names = (
    "organism_taxon" => 'Genome::Taxon',
    "organism_individual" => 'Genome::Individual',
    "population_group" => 'Genome::PopulationGroup',
    "organism_sample" => 'Genome::Sample',
    "sample_attribute" => 'Genome::SubjectAttribute',
    "population_group_member" => 'Genome::SubjectAttribute',
);
sub genome_class_name {
    my $self = shift;

    my $lims_table_name = $self->lims_table_name;
    return if not $lims_table_name;

    my $genome_class_name = $subject_class_names_to_genome_class_names{$lims_table_name};
    if ( not $genome_class_name ) {
        $self->error_message('No genome class name for lims table name => '.$lims_table_name);
        return;
    }

    return $genome_class_name;
}

sub genome_entity {
    my $self = shift;

    my $genome_class_name = $self->genome_class_name;
    return if not $genome_class_name;

    my $subject_id = $self->subject_id;
    if($genome_class_name eq 'Genome::SubjectAttribute') {
        $self->error_message("Cannot use single misc update to get genome entity for subject class name (".$self->subject_class_name."). Use MultiMiscUpdate.");
        return;
    }

    my $genome_entity = $genome_class_name->get($subject_id);
    if ( not $genome_entity ) {
        # Should exist
        $self->error_message("Failed to get $genome_class_name for id => $subject_id");
        return;
    }

    return $genome_entity;
}

my %subject_class_names_to_site_tgi_class_names = (
    "organism_taxon" => 'Genome::Site::TGI::Synchronize::Classes::Taxon',
    "organism_individual" => 'Genome::Site::TGI::Synchronize::Classes::Individual',
    "population_group" => 'Genome::Site::TGI::Synchronize::Classes::PopulationGroup',
    "organism_sample" => 'Genome::Site::TGI::Synchronize::Classes::Sample',
    "sample_attribute" => 'Genome::Site::TGI::Synchronize::Classes::SubjectAttribute',
    "population_group_member" => 'Genome::Site::TGI::Synchronize::Classes::PopulationGroupMember',
);
sub perform_update {
    my $self = shift;

    my $lims_table_name = $self->lims_table_name;
    return if not $lims_table_name;

    my $site_tgi_class_name = $subject_class_names_to_site_tgi_class_names{$lims_table_name};
    if ( not $site_tgi_class_name ) {
        $self->error_message('No site tgi class name for lims table name => '.$lims_table_name);
        return $self->failure;
    }

    my $lims_property_name = $self->subject_property_name;
    my $genome_property_name = $site_tgi_class_name->lims_property_name_to_genome_property_name($lims_property_name);
    if ( not $genome_property_name ) {
        $self->error_message('No genome property name for lims property name => '.$lims_property_name);
        return $self->failure;
    }

    if ( grep { $lims_table_name eq $_ } (qw/ sample_attribute population_group_member /) ) {
        $self->error_message("Cannot UPDATE $lims_table_name! It must be INSERTED then DELETED!");
        return $self->failure;
    }

    if ( not grep { $genome_property_name eq $_ } $site_tgi_class_name->properties_to_keep_updated ) {
        $self->status_message('Update for genome property name not supported => '.$genome_property_name);
        return $self->skip;
    }

    my $genome_class_name = $self->genome_class_name;
    if ( not $genome_class_name ) {
        return $self->failure; 
    }

    my $genome_entity = $self->genome_entity;
    if ( not $genome_entity ) {
        return $self->failure;
    }

    # Get values
    my $new_value = $self->new_value;
    my $old_value = $self->old_value;
    my ($current_attr, $current_value);
    if ( $genome_property_name eq 'name' ) {
        $current_value = $genome_entity->name;
    }
    else {
        $current_attr = Genome::SubjectAttribute->get(
            subject_id => $genome_entity->id,
            attribute_label => $genome_property_name,
            nomenclature => 'WUGC',
        );
        $current_value = $current_attr->attribute_value if $current_attr;
    }
    $self->{_current_value} = $current_value;

    # NEW and CURRENT are NULL
    if ( not defined $new_value and not defined $current_value ) {
        $self->is_reconciled(1);
        return $self->success;
    }

    # NEW and CURRENT are NOT NULL and the same
    if ( defined $current_value and defined $new_value and $current_value eq $new_value ) {
        $self->is_reconciled(1);
        return $self->success;
    }

    # OLD and CURRENT value do not match
    if ( defined $old_value and defined $current_value and $old_value ne $current_value ) {
        $self->error_message("Current APipe value ($current_value) does not match the LIMS old value ($old_value)!");
        return $self->failure;
    }

    # Update
    my $updated_value;
    if ( $genome_property_name eq 'name' ) {
        $updated_value = $genome_entity->name($new_value);
    }
    else {
        $current_attr->delete if $current_attr;
        my $new_attr = Genome::SubjectAttribute->create(
            subject_id => $genome_entity->id,
            attribute_label => $genome_property_name,
            attribute_value => $new_value,
            nomenclature => 'WUGC',
        );
        $updated_value = $new_attr->attribute_value
    }

    if ( not $updated_value or $updated_value ne $new_value ) {
        $self->error_message('Failed to set new value!');
        return $self->failure;
    }

    return $self->success;
}

sub _set_result {
    my ($self, $result) = @_;

    $self->result($result);
    $self->status(
        join("\t", 
            $result, 
            map({ $self->$_; } (qw/ subject_class_name subject_id subject_property_name /)),
            "'".( $self->{_current_value} // 'NA' )."'",
            map({ defined $_ ? "'".$_."'" : "'NULL'"; }  $self->old_value, $self->new_value,),
        )
    );
}

sub success {
    my $self = shift;
    $self->_set_result( uc($self->description) );
    $self->is_reconciled(1);
    return 1;
}

sub skip {
    my $self = shift;
    $self->_set_result('SKIP');
    return 1;
}

sub failure {
    my $self = shift;
    $self->_set_result('FAILED');
    return;
}

1;

