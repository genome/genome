package Genome::Subject;

use strict;
use warnings;
use Genome;

class Genome::Subject {
    is => 'Genome::Notable',
    table_name => 'subject.subject',
    is_abstract => 1,
    subclassify_by => 'subclass_name',
    id_generator => '-uuid',
    id_by => [
        subject_id => {
            is => 'Number',
            len => 10,
        },
    ],
    has => [
        subclass_name => {
            is => 'Text',
            len => 255,
        },
        subject_type => {
            is_abstract => 1,
            is_transient => 1,
            doc => 'Plain text description of the type of subject, defined in subclasses',
        },
    ],
    has_optional => [
        name => {
            is => 'Text',
            len => 255,
            doc => 'Official name of the subject',
        },
        common_name => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            is_mutable => 1,
            is_many => 0,
            where => [ attribute_label => 'common_name' ],
        },
        description => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            is_mutable => 1,
            is_many => 0,
            where => [ attribute_label => 'description' ],
            doc => 'Some plain-text description of the subject',
        },
        nomenclature => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            is_mutable => 1,
            is_many => 0,
            where => [ attribute_label => 'nomenclature', nomenclature => 'WUGC' ],
            doc => 'The nomenclature that information about this subject follows (TCGA, WUGC, etc)',
        },
    ],
    has_many_optional => [
        attributes => {
            is => 'Genome::SubjectAttribute',
            reverse_as => 'subject',
        },
        project_parts => {
            is => 'Genome::ProjectPart',
            reverse_as => 'entity',
            is_mutable => 1,
        },
        projects => {
            is => 'Genome::Project',
            via => 'project_parts',
            to => 'project',
            is_mutable => 1,
            doc => 'Projects that include this subject.',
        },
        project_names => {
            is => 'Text',
            via => 'projects',
            to => 'name',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'Contains all information about a particular subject (sample, individual, etc)',
};

# TODO Figure out a better name
# This method should be overridden in subclasses and should return the object "further up"
# the derivation tree. For example, a sample would override this to return the individual from
# which it is derived.
sub get_source {
    return;
}

sub get_source_with_class {
    my ($self, $class) = @_;
    my $source = $self;
    while ($source) {
        $source = $source->get_source;
        last unless $source;
        return $source if $source->isa($class);
    }
    return;
}

sub attributes_for_nomenclature {
    my ($self, $nom) = @_;
    die "must supply nomenclature" if !$nom;
    my @fields = $nom->fields();
    my @attr = Genome::SubjectAttribute->get(
        nomenclature => [ map {$_->id} @fields ],
        subject_id => $self->id
    );
    return @attr;
}

sub __display_name__ {
    my $self = shift;
    my $name = $self->name . ' ' if defined $self->name;
    my $common_name = ($self->can("common_name") ? $self->common_name : ());
    $name .= '(' . ($common_name ? $common_name . ' ' : '') . $self->id . ')';
    return $name;
}

sub create {
    my $class = shift;

    my ($bx, %extra) = $class->define_boolexpr(@_);
    return if not $bx;

    my $self = $class->SUPER::create($bx);
    return if not $self;

    my $nomenclature = $self->nomenclature;
    $nomenclature = 'WUGC' if not defined $nomenclature;
    for my $label (sort keys %extra) {
        my $attribute = Genome::SubjectAttribute->create(
            attribute_label => $label,
            attribute_value => $extra{$label},
            subject_id => $self->id,
            nomenclature => $nomenclature,
        );
        unless ($attribute) {
            $self->error_message("Could not create attribute $label => ".$extra{$label}."!");
            $self->delete;
            return;
        }
    }

    return $self;
}

sub delete {
    my $self = shift;
    my @attributes = $self->attributes;
    for my $attribute (@attributes) {
        Carp::confess "Could not delete attribute " . $attribute->attribute_label . " for subject " . $self->id unless $attribute->delete;
    }
    return $self->SUPER::delete;
}

1;

