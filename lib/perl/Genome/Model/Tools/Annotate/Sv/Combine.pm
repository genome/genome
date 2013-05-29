package Genome::Model::Tools::Annotate::Sv::Combine;

use strict;
use warnings;
use Genome;

my @has_param;
my %instances;

BEGIN{
    my @annotators = (
        'Genome::Model::Tools::Annotate::Sv::Transcripts',
        'Genome::Model::Tools::Annotate::Sv::FusionTranscripts',
        'Genome::Model::Tools::Annotate::Sv::Dbsnp',
        'Genome::Model::Tools::Annotate::Sv::Segdup',
        'Genome::Model::Tools::Annotate::Sv::RepeatMasker',
        'Genome::Model::Tools::Annotate::Sv::Dbvar',
    );
    foreach my $module (@annotators) {
        my $module_meta = UR::Object::Type->get($module);
        my @module_path = split /::/, $module;
        my $module_short_name = Genome::Utility::Text::camel_case_to_string($module_path[-1]);
        $module_short_name =~ s/\s+/_/g;
        my @p = $module_meta->properties;
        foreach my $p (@p) {
            if ($p->can("is_input") and $p->is_input) {
                my $name = $p->property_name;
                $name = $module_short_name."_".$name;
                my %data = %{$p};
                for my $key (keys %data) {
                    delete $data{$key} if $key =~ /^_/;
                }
                delete $data{id};
                delete $data{db_committed};
                delete $data{class_name};
                delete $data{is_input};
                $data{is_optional} = 1;
                $data{is_param} = 1;
                $data{doc} .= " -- for the $module_short_name annotator";
                $data{property_name} = $name;
                push @has_param, $name, \%data;
            }
        }
    }
}

class Genome::Model::Tools::Annotate::Sv::Combine {
    is => "Genome::Model::Tools::Annotate::Sv::Base",
    doc => "Combine multiple SV annotation methods",
    has_param => \@has_param,
    has_input => [
        annotator_list  => {
        is => 'String',
        is_many => 1,
        default => ['Transcripts', 'FusionTranscripts', 'Dbsnp'],
        },
    ],
};

sub help_detail {
    return "Run one or more SV annotators and combine their output.  Specify which
    annotators to run with --annotator-list (e.g. --annotator-list=Transcripts,FusionTranscripts,Dbsnp).
    Parameters for the individual annotators are prefixed with the annotator name (i.e. the --fusion-output-file
    parameter for the FusionTranscript annotator should be specified with --fusion-transcripts-fusion-output-file";
}

sub process_breakpoint_list {
    my $self = shift;
    my $breakpoints_list = shift;
    my %all_content;

    for my $type ($self->annotator_list) {
        unless ($instances{$type}) {
            $instances{$type} = $self->_create_instance_of_type($type);
        }
        my $content = $instances{$type}->process_breakpoint_list($breakpoints_list);
        foreach my $key (keys %$content) {
            push @{$all_content{$key}}, @{$content->{$key}};
        }
    }
    return \%all_content;
}

sub column_names {
    my $self = shift;
    my @all_column_names;
    foreach my $type ($self->annotator_list) {
        unless ($instances{$type}) {
            $instances{$type} = $self->_create_instance_of_type($type);
        }
        @all_column_names = (@all_column_names, $instances{$type}->column_names);
    }
    return @all_column_names;
}

sub _create_instance_of_type {
    my $self = shift;
    my $type = shift;
    my $class_name = "Genome::Model::Tools::Annotate::Sv::$type";
    my $module_meta = UR::Object::Type->get($class_name);
    my $module_short_name = Genome::Utility::Text::camel_case_to_string($type);
    $module_short_name =~ s/\s+/_/g;
    my @p = $module_meta->properties;
    my %params;
    foreach my $property (@p) {
        if ($property->can("is_input") and $property->is_input) {
            my $property_name = $module_short_name."_".$property->property_name;
            if ($self->$property_name) {
                $params{$property->property_name} = $self->$property_name;
            }
        }
    }
    $params{input_file} = $self->input_file;
    $params{output_file} = $self->output_file;
    $params{annotation_build} = $self->annotation_build;
    my $instance = $class_name->create(%params);
    return $instance;
}


1;

