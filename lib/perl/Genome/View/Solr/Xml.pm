package Genome::View::Solr::Xml;

use strict;
use warnings;
use Genome;

use WebService::Solr;
use MIME::Base64;

class Genome::View::Solr::Xml {
    is => 'UR::Object::View::Default::Xml',
    is_abstract => 1,
    attributes_have => [
        is_field => {
            is => 'Boolean',
            default_value => 0,
        },
        field_name => {
            is => 'Text',
            default_value => '',
        },
    ],
    has => [
        _doc    => {
            is_transient => 1,
            doc => 'the WebService::Solr::Document object used to generate this view'
        },
    ],
    has_field => [
        type => {
            is_abstract => 1,
            is => 'Text',
            doc => 'The type represented by the document as referred to in Solr--override this in subclasses'
        },
        subject_class => {
            is => 'Text',
            via => 'subject',
            to => 'class',
            field_name => 'class',
        },
        old_solr_id => {
            is => 'Text',
            field_name => 'id',
            calculate => q| $self->_generate_id_field_data |,
        },
        title => {
            is => 'Text',
            calculate => q| $self->_generate_title_field_data |,
        },
        timestamp => {
            is => 'Text',
            calculate => q| $self->_generate_timestamp_field_data |,
        },
        object_id => {
            is => 'Text',
            calculate => q| $self->_generate_object_id_field_data |,
        },
    ],
    has_constant => [
        perspective => {
            value => 'solr',
        },
    ],
    doc => 'The base class for creating the XML document that Solr indexes for an object.'
};

sub create {
    my $class = shift;
    my %params = @_;

    if(exists $params{content_doc}) {
        return $class->_reconstitute_from_doc($params{content_doc});
    } else {
        return $class->SUPER::create(%params);
    }
}

sub _reconstitute_from_doc {
    my $class = shift;
    my $solr_doc = shift;

    unless($solr_doc->isa('WebService::Solr::Document')) {
        $class->error_message('content_doc must be a WebService::Solr::Document');
        return;
    }

    my $subject = Genome::Search->get_subject_from_doc($solr_doc);
    unless ($subject) {
        die 'Failed to get subject from solr_doc.';
    }
    my $self = $class->SUPER::create(subject_id => $subject->id, subject_class_name => $subject->class);

    $self->_doc($solr_doc);

    my $widget = $self->widget();
    my ($content_ref,$fh) = @$widget;
    $$content_ref = $self->_doc->to_xml;

    return $self;
}


sub content_doc {
    my $self = shift;
    my $content = $self->content; #force document generation
    return $self->_doc;
}

sub _generate_fields {
    my $self = shift;
    my @fields;
    my $super_content = $self->SUPER::_generate_content;
    my @text_field_properties = $self->__meta__->properties(is_field => 1, data_type => 'Text');
    for my $text_field_property (@text_field_properties) {
        my $key = ($text_field_property->field_name || $text_field_property->property_name);
        my $property_name = $text_field_property->property_name;
        my $value = $self->$property_name || '';
        push @fields, [$key, $value];
    }
    my @aspect_sets = $self->aspects(-group_by => 'position');
    for my $aspect_set (@aspect_sets) {
        my @aspects = $aspect_set->members;
        my $key = $aspects[0]->position;
        if (grep { $key eq $_ } ('timestamp', 'title')) { # these are calculated and will have already inspected the aspect
            next;
        }
        my @values;
        for my $aspect (@aspects) {
            if ($aspect->delegate_view && $aspect->delegate_view->isa('Genome::View::Solr::Xml')) {
                my $prefix = $aspect->name;
                my @child_fields = $aspect->delegate_view->_generate_fields;
                for my $child_field (@child_fields) {
                    my ($child_key, $child_value) = @$child_field;
                    push @fields, [$child_key, $child_value, $prefix];
                }
            } else {
                my $aspect_content = $self->_generate_content_for_aspect($aspect);
                my $node_list = $aspect_content->findnodes('//value');
                if ($node_list->isa('XML::LibXML::NodeList')) {
                    if ($key eq 'content') {
                        push @values, $aspect->name . ':';
                    }
                    while (my $node = $node_list->shift) {
                        push @values, $node->to_literal;
                    }
                } else {
                    my $value = $aspect_content->findvalue('value');
                    if (not defined $value) {
                        die "$key has an undefined value and no delegate view";
                    }
                    push @values, $value;
                }
            }
        }
        my $value = join(' ', @values);
        push @fields, [$key, $value];
    }
    return @fields;
}

sub _generate_content {
    my $self = shift;

    my $subject = $self->subject;
    return unless $subject;

    my @fields = $self->_generate_fields;

    my @schema_fields = qw(
        id
        idAnalyzed
        object_id
        class
        type
        title
        content
        display_title
        display_type
        display_icon_url
        display_url0
        display_label1
        display_url1
        display_label2
        display_url2
        display_label3
        display_url3
        word
        timestamp
    );
    my @excluded_content_keys = grep { $_ ne 'content' } @schema_fields;
    my @content_values;
    my @solr_fields;
    for my $field (@fields) {
        my ($key, $value, $prefix) = @$field;
        my $prefixed_key = ($prefix ? join('_', $prefix, $key) : $key);

        unless (grep { $prefixed_key eq $_ } @schema_fields) {
            unless (grep { $key =~ /^$_$/ } @excluded_content_keys) {
                push @content_values, "$key:$value"; # dynamic fields are not yet searchable so pump into content field
            }
            $prefixed_key .= "_t"; # required for dynamic fields. later should add solr_type attribute to map to non-text
        }

        if ($key eq 'content') { # building up a single content field
            push @content_values, $value;
        } else {
            push @solr_fields, WebService::Solr::Field->new($prefixed_key => $value);
        }
    }
    push @solr_fields, WebService::Solr::Field->new(content => join(' ', @content_values));

    $self->_doc( WebService::Solr::Document->new(@solr_fields) );
    return $self->_doc->to_xml;
}

sub _generate_title_field_data {
    my $self = shift;
    my $subject = $self->subject;

    my @aspects = $self->aspects;
    my @title_aspects = grep($_->position eq 'title', @aspects);

    my $title;

    if(scalar @title_aspects) {
        my @title_parts;

        for my $aspect (@title_aspects) {
            my $property = $aspect->name;

            my $value = $subject->$property;

            $value = '' unless defined $value; #Not useful for indexing, but will add an extra space in case anyone cares.

            push @title_parts, $value;
        }

        $title = join(' ', @title_parts);
    }

    unless($title) {
        if($subject->can('name') and $subject->name) {
            $title = $subject->name;
        } else {
            $title = $self->type . ' ' . $subject->id;
        }
    }

    return $title;
}

sub _generate_id_field_data {
    my $self = shift;
    my $subject = $self->subject;

    #TODO after all code is setting solr_id, make that field the primary id for Solr;
    #then remove the class and '---' from this field
    return $subject->class . '---' . $self->_generate_object_id_field_data;
}

sub _generate_object_id_field_data {
    my $self = shift;
    my $subject = $self->subject;
    my $object_id = $subject->id;

    # Sets have invalid XML chars in their IDs so we encode them in Base64.
    # Decoding is done in Genome::Search::get_subject_from_doc so keep symmetry there.
    # TODO Encoding/decoding should probably be handled by the object itself.
    if ($subject->isa('UR::Object::Set')) {
        $object_id = encode_base64($object_id);
    }

    return $object_id;
}

sub _generate_timestamp_field_data {
    my $self = shift;
    my $subject = $self->subject;

    my @aspects = $self->aspects;
    my @timestamp_aspects = grep($_->position eq 'timestamp', @aspects);

    #By convention this timestamp is used when we don't know the real timestamp
    my $timestamp = '1999-01-01T01:01:01Z';

    if(scalar @timestamp_aspects) {
        if(scalar @timestamp_aspects > 1) {
            $self->error_message('Only one timestamp may be supplied.');
            die $self->error_message;
        }
        my $aspect = $timestamp_aspects[0];
        my $property = $aspect->name;

        my $value = $subject->$property;

        if($value) {
            my $solr_timestamp_format = qr(\d{4}-\d{1,2}-\d{1,2}T\d{1,2}:\d{1,2}:\d{1,2}Z);

            if($value =~ $solr_timestamp_format) {
                $timestamp = $value;
            } else {
                my ($a, $b) = split / /, $value;
                $b =~ s/.\d{6}$//g;
                $timestamp = sprintf("%sT%sZ", $a, $b);
            }

            unless($timestamp =~ '\d{4}-\d{1,2}-\d{1,2}T\d{1,2}:\d{1,2}:\d{1,2}Z') {
                $self->error_message('Could not parse timestamp: ' . $timestamp . '(format should be yyyy-mm-ddThh:mm:ssZ)');
                die $self->error_message;
            }
        }
    }

    return $timestamp;
}

1;
