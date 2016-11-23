package Genome::InstrumentData::Command::List;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::List {
    is => 'Genome::Object::Command::List',
    has => [
        subtype => {
            is => 'Text',
            default_value => 'any',
            is_optional => 1,
            shell_args_position => 1,
            valid_values => [qw(454 sanger solexa imported any)],
            doc => 'Filter by subclass of InstrumentData.  Each provides different default values of the "show" parameter.',
        },

        filter => {
            is => 'Text',
            shell_args_position => 2,
            doc => 'Filter results based on the parameters.  See below for details.',
        },
    ],
    doc => 'find instrument data',
};

my %CLASS_LOOKUP = (
    '454' => 'Genome::InstrumentData::454',
    imported => 'Genome::InstrumentData::Imported',
    solexa => 'Genome::InstrumentData::Solexa',
    any => 'Genome::InstrumentData',
);
sub subject_class_name {
    my $self = shift;

    return $CLASS_LOOKUP{$self->subtype};
}

sub _resolve_field_list {
    my $self = shift;

    my $default_show = $self->_default_show;
    $self->_set_show_default_value($default_show);

    unless(defined $self->show) {
        $self->show($default_show);
    }

    return $self->SUPER::_resolve_field_list(@_);
}

sub _set_show_default_value {
    my $self = shift;
    my $show_value = shift;

    $self->__meta__->property('show')->default_value($show_value);
}


my %DEFAULT_SHOW = (
    '454' => [qw(
        id
        run_name
        region_number
        index_sequence
        sample_name
    )],

    solexa => [qw(
        id
        flow_cell_id
        lane
        index_sequence
        sample_name
        library_name
        read_length
        is_paired_end
        clusters
        target_region_set_name
    )],

    imported => [qw(
        id
        sample_name
        sequencing_platform
        import_format
    )],

    any => [qw(
        id
        sample_name
        library_name
        run_name
        subset_name
    )],
);
sub _default_show {
    my $self = shift;

    return join(',', @{$DEFAULT_SHOW{$self->subtype}});
}

my %BASE_FILTER_LOOKUP = (
    sanger => 'sequencing_platform=sanger',
    solexa => 'sequencing_platform=solexa',
    '454' => 'sequencing_platform=454',
);
sub _base_filter {
    my $self = shift;

    if (exists $BASE_FILTER_LOOKUP{$self->subtype}) {
        return $BASE_FILTER_LOOKUP{$self->subtype};
    } else {
        return;
    }
}

sub _hint_string {
    my $self = shift;

    my @show_parts = grep { $self->_show_item_is_property_name($_) } $self->_resolve_field_list();
    my $class_name = $self->subject_class_name;

    my @hint_parts = ();
    my $include_attributes = 0;

    for my $part (@show_parts) {
        my $prop = $class_name->__meta__->property(property_name => $part);
        if($prop and (($prop->via and $prop->via eq 'attributes') or $prop->is_calculated)) {
            $include_attributes = 1;
        } else {
            push @hint_parts, $part;
        }
    }
    push @hint_parts, 'attributes' if $include_attributes;

    return join(',', @hint_parts);

}

1;
