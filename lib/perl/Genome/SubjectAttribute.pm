package Genome::SubjectAttribute;

use strict;
use warnings;
use Genome;

class Genome::SubjectAttribute {
    table_name => 'subject.subject_attribute',
    id_by => [
        attribute_label => {
            is => 'Text',
            len => 64,
            default_value => 'NONE',
        },
        attribute_value => {
            is => 'Text',
            len => 512,
        },
        subject_id => {
            is => 'Number',
            len => 10,
        },
        nomenclature => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        nomenclature_field => {
            is => 'Genome::Nomenclature::Field',
            id_by => 'nomenclature',
        },
        subject => {
            is => 'Genome::Subject',
            id_by => 'subject_id',
        },
        subject_name => {
            via => 'subject',
            to => 'name',
        },
        nomenclature_field_name => {
            via => 'nomenclature_field',
            to => 'name',
        },
        nomenclature_field_type => {
            via => 'nomenclature_field',
            to => 'type',
        },
        nomenclature_name => { via => 'nomenclature_field' },
        nomenclature_id => { via => 'nomenclature_field' },
        nomenclature_obj => {
            is => 'Genome::Nomenclature',
            via => 'nomenclature_field',
            to => 'nomenclature',
        },
        all_nomenclature_fields => {
            via => 'nomenclature_obj',
            to => 'fields',
        },
    ],
    has_optional => [
        _individual => {
            is => 'Genome::Individual',
            id_by => 'attribute_value',
        },
        sample => {
            is => 'Genome::Sample',
            id_by => 'subject_id',
            constraint_name => 'GSA_GS_FKPK',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'Represents a particular attribute of a subject',
};

sub create {
    my $class = shift;
    my $bx = $class->define_boolexpr(@_);
    # TODO This is a workaround that allows nomenclature to be in the id_by block
    # and have a default value. Doing so in the class definition doesn't work due
    # to some sort of UR bug that Tony is aware of.
    unless ($bx->specifies_value_for('nomenclature')) {
        $bx = $bx->add_filter('nomenclature' => 'WUGC');
    }
    return $class->SUPER::create($bx);
}


1;

