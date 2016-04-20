package Genome::ModelGroup::View::Status::Xml;

use strict;
use warnings;
use Genome;
use Data::Dumper;
use XML::LibXML;

class Genome::ModelGroup::View::Status::Xml {
    is => 'UR::Object::View::Default::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'id',
                'name',
                'user_name',
                'model_count',
                {
                    name => 'models',
                    subject_class_name => 'Genome::Model',
                    aspects => [
                        'name',
                        'run_as',
                        'created_by',
                        'creation_date',
                        {
                            name => 'subject',
                            aspect => [
                                'id', 'name', 'subclass_name'
                            ],
                            perspective => 'default',
                            toolkit => 'xml',
                            subject_class_name => 'Genome::Subject'
                        },
                        {
                            name => 'processing_profile',
                            aspects => [
                                'id', 'name', 'type_name', 'subclass_name',
                            ],
                            perspective => 'default',
                            toolkit => 'xml',
                            subject_class_name => 'Genome::ProcessingProfile',
                        },
                        {
                            name => 'inputs',
                            aspects => [
                                'value_id',
                                'name',
                                'value',
                                'value_class_name',
                            ],
                            perspective => 'default',
                            toolkit => 'xml',
                        },
                        {
                            name => 'builds',
                            aspects => [
                                'id',
                                'data_directory',
                                'status',
                                'date_scheduled',
                                'date_completed',
                                {
                                    name => 'notes',
                                    aspects => [
                                        'editor_id',
                                    ],
                                    perspective => 'default',
                                    toolkit => 'xml',
                                },
                                {
                                    name => 'delta_model_input_differences_from_model',
                                    aspects => ['value_id', 'name', 'value', 'value_class_name'],
                                    perspective => 'default',
                                    toolkit => 'xml',
                                    subject_class_name => 'Genome::Model::Input',
                                },
                                {
                                    name => 'build_input_differences_from_model',
                                    aspects => ['value_id', 'name', 'value', 'value_class_name'],
                                    perspective => 'default',
                                    toolkit => 'xml',
                                    subject_class_name => 'Genome::Model::Build::Input',
                                },
                            ],
                            perspective => 'default',
                            toolkit => 'xml',
                            subject_class_name => 'Genome::Model::Build',
                        },
                    ],
                    perspective => 'default',
                    toolkit => 'xml',
                },
            ]
        }
    ]
};

sub _generate_content {
    my $self = shift;
    my $subject = $self->subject;

    if($subject) {
        #preload data for efficiency
        my @m = $subject->models;
        my @in = Genome::Model::Input->get(model_id => [map($_->id, @m)]);
        my @i = Genome::InstrumentData->get([map($_->value_id, grep { $_->name eq 'instrument_data' } @in)]);
        Genome::InstrumentDataAttribute->get(instrument_data_id => [map($_->id, @i)], attribute_label => 'flow_cell_id');
        Genome::Model::Link->get(from_model_id => [map($_->id, @m)]);
        Genome::Model::Link->get(to_model_id => [map($_->id, @m)]);
        my @b = Genome::Model::Build->get(model_id => [map($_->id, @m)]);
        Genome::Model::Build::Input->get(build_id => [map($_->id, @b)]);
        Genome::Model::Event->get(build_id => [map($_->id, @b)]);
        Genome::MiscNote->get(subject_id => [map($_->id, (@m, @b))]);
        my @s = Genome::Subject->get(id => [map($_->subject_id, @m)]);
        Genome::SubjectAttribute->get(subject_id => [map($_->id, @s)]);
    }

    return $self->SUPER::_generate_content(@_);
}

1;

=pod

=head1 NAME

Genome::ModelGroup::View::Status::XML - status summary for a model group in XML format

=head1 SYNOPSIS

$i = Genome::ModelGroup->get(1234);
$v = Genome::ModelGroup::View::Status::Xml->create(subject => $i);
$xml = $v->content;

=head1 DESCRIPTION

This view renders the summary of an model group's status in XML format.

=cut

