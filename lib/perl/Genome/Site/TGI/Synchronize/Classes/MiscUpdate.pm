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
    table_name => 'misc_update',
    is => 'UR::Object',
    is_abstract => 1,
    subclassify_by => 'subclass_name',
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
            is => 'Text', 
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
    has_optional_calculated_constant => [
        schema => {
            calculate_from => 'subject_class_name',
            calculate => sub {
                my $subject_class_name = shift;
                return (schema_and_lims_table_name($subject_class_name))[0];
            },
        },
        lims_table_name => {
            calculate_from => 'schema_and_lims_table_name',
            calculate => sub {
                my $subject_class_name = shift;
                return (schema_and_lims_table_name($subject_class_name))[1];
            },
        },
        subclass_name => {
            is => 'Text',
            calculate_from => 'subject_class_name',
            calculate => sub {
                my $subject_class_name = shift;
                my $lims_table_name = schema_and_lims_table_name_from_subject_class_name($subject_class_name);
                my $class = __PACKAGE__.'::'.Genome::Utility::Text::string_to_camel_case($lims_table_name);
                return $class if eval{ $class->__meta__; };
                return __PACKAGE__.'::Unsupported';
            },
        },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

sub schema_and_lims_table_name_from_subject_class_name {
    my $subject_class_name = shift;
    Carp::confess('No subject class name to get schema and lims table name') if not $subject_class_name;
    my ($schema, $lims_table_name) = split(/\./, $subject_class_name);
    if ( not $schema ) {
        Carp::confess('Failed to get schema from subject class name => '.$subject_class_name);
    }
    if ( not $lims_table_name ) {
        Carp::confess('Failed to get lims table name from subject class name => '.$subject_class_name);
    }
    return ( $schema, $lims_table_name );
}

my %lims_table_names_to_site_tgi_class_names = (
    "organism_taxon" => 'Genome::Site::TGI::Synchronize::Classes::OrganismTaxon',
    "organism_individual" => 'Genome::Site::TGI::Synchronize::Classes::OrganismIndividual',
    "population_group" => 'Genome::Site::TGI::Synchronize::Classes::PopulationGroup',
    "organism_sample" => 'Genome::Site::TGI::Synchronize::Classes::OrganismSample',
    "library_summary" => 'Genome::Site::TGI::Synchronize::Classes::LibrarySummary',
    "index_illumina" => 'Genome::Site::TGI::Synchronize::Classes::IndexIllumina',
    "region_index_454" => 'Genome::Site::TGI::Synchronize::Classes::RegionIndex454',
    "run_region_454" => 'Genome::Site::TGI::Synchronize::Classes::RegionIndex454',
);
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

sub site_tgi_class_name {
    my $self = shift;

    my $lims_table_name = $self->lims_table_name;
    my $site_tgi_class_name = $lims_table_names_to_site_tgi_class_names{$lims_table_name};
    if ( not $site_tgi_class_name ) {
        $self->error_message('Unsupported LIMS table name => '.$lims_table_name);
        return;
    }

    return $site_tgi_class_name;
}

my %subject_class_names_to_genome_class_names = (
    "organism_taxon" => 'Genome::Taxon',
    "organism_individual" => 'Genome::Individual',
    "population_group" => 'Genome::PopulationGroup',
    "organism_sample" => 'Genome::Sample',
    "sample_attribute" => 'Genome::SubjectAttribute',
    "population_group_member" => 'Genome::SubjectAttribute',
    "library_summary" => 'Genome::Library',
    "index_illumina" => 'Genome::InstrumentData::Solexa',
    "region_index_454" => 'Genome::InstrumentData::454',
    "run_region_454" => 'Genome::InstrumentData::RunRegion454',
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

    return $self->{_genome_entity} if $self->{_genome_entity};

    my $genome_class_name = $self->genome_class_name;
    return $self->failure if not $genome_class_name;

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
    $self->{_genome_entity} = $genome_entity;

    return $genome_entity;
}

sub genome_property_name {
    my $self = shift;

    my $site_tgi_class_name = $self->site_tgi_class_name;
    return if not $site_tgi_class_name;

    my $subject_property_name = lc( $self->subject_property_name );
    my $genome_property_name = $site_tgi_class_name->lims_property_name_to_genome_property_name($subject_property_name);
    if ( not $genome_property_name ) {
        $self->error_message("Failed to get genome property name! $site_tgi_class_name => ".$self->subject_property_name);
        return;
    }

    return $genome_property_name;
}

sub _set_result {
    my ($self, $result) = @_;

    $self->result($result);
    $self->status(
        join("\t", 
            $result, 
            map({ $self->$_; } (qw/ description subject_class_name subject_id /)),
            lc( $self->subject_property_name ),
            "'".( $self->{_current_value} // 'NA' )."'",
            map({ defined $_ ? "'".$_."'" : "'NULL'"; }  $self->old_value, $self->new_value,),
        )
    );
}

sub success {
    my $self = shift;
    $self->_set_result('PASS');
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
    $self->_set_result('FAIL');
    return;
}

sub has_failed {
    my $self = shift;
    return 1 if $self->result and $self->result eq 'FAIL';
}

1;

