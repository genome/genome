package Genome::Model::SomaticValidation::Command::ImportVariants;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::ImportVariants {
    is => 'Command::V2',
    has_input => [
        variant_file_list => {
            is => 'Text',
            doc => 'File listing the variants to be uploaded',
        },
        models => {
            is => 'Genome::Model',
            doc => 'The somatic variation models for the variants in the list, e.g., as a group (model_groups.id=?) or comma-delimited list of ids',
            is_many => 1,
        },
        variant_file_format => {
            is => 'Text',
            doc => 'format of the files in the variant_file_list',
            default_value => 'annotation',
            is_optional => 1,
        },
    ],
    has_optional_output => [
        results => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            is_many => 1,
            doc => 'The results from running this command',
        },
,
    ],
    doc => 'enter the variants for validation from many previous somatic-variation runs at once',
};

sub sub_command_category { 'analyst tools' }

sub execute {
    my $self = shift;

    #First, find all the files we'll be working with
    my %data;

    my $variant_file_list_fh = Genome::Sys->open_file_for_reading($self->variant_file_list);

    my $variant_type;
    while(my $line = <$variant_file_list_fh>) {
        chomp $line;
        if($line =~ m{.*/((?:\w+\d+)|(?:\w_\w\w-[^/]+))/[^/]+$}) {
            my $patient = $1;

            $data{$patient}{$variant_type} ||= [];
            push @{ $data{$patient}{$variant_type} }, $line;
        } elsif($line =~ m{.*/([^/.]+)(?:\.[^./]+)+$}) {
            my $patient = $1;

            $data{$patient}{$variant_type} ||= [];
            push @{ $data{$patient}{$variant_type} }, $line;
        } elsif(grep($_ eq $line, 'snvs', 'indels', 'svs')) {
            $variant_type = $line; #header indicating SNVs, indels, or SVs
            $variant_type =~ s/s$//;
        } else {
            die $self->error_message('Could not determine patient for this file: ' . $line);
        }
    }

    my @results;
    for my $m ($self->models) {
        my $patient = $m->subject;
        if($patient->isa('Genome::Sample')) { $patient = $patient->source; }
        unless($patient) {
            die $self->error_message('No patient found linked to subject of model ' . $m->__display_name__);
        }

        my $model_found = 0;
        ACCESSOR: for my $accessor ('id', 'name', 'common_name') {
            if(exists $data{$patient->$accessor}) {
                for my $variant_type (keys %{$data{$patient->$accessor}}) {
                    my $result = $self->_upload_result($m, $variant_type, @{ delete $data{$patient->$accessor}{$variant_type} });
                    push @results, $result;
                }
                $model_found = 1;
                last ACCESSOR;
            }
        }

        unless($model_found) {
            $self->warning_message('No data based on model: ' . $m->__display_name__);
        }
    }

    for my $patient (keys %data) {
        $self->warning_message('No model found to create results for ' . $patient . '.');
    }

    $self->results(\@results);
    $self->status_message('Successfully created results.');
    $self->status_message('Result IDs: ' . join(',', map($_->id, @results)));

    return 1;
}

sub _upload_result {
    my $self = shift;
    my $model = shift;
    my $variant_type = shift;
    my @files = @_;

    my $file;
    if(scalar(@files) > 1) {
        $file = Genome::Sys->create_temp_file_path;
        Genome::Sys->cat(
            input_files => \@files,
            output_file => $file,
        );
    } elsif (scalar(@files) == 1) {
        $file = $files[0];
    } else {
        die '_upload_result called with no variant files';
    }

    my $build = $model->last_complete_build;
    unless($build) {
        die $self->error_message('No complete build found for model ' . $model->__display_name__);
    }

    my $result_cmd = Genome::Model::SomaticValidation::Command::ManualResult->create(
        source_build => $build,
        variant_file => $file,
        variant_type => $variant_type,
        format => $self->variant_file_format,
        description => 'generated from ' . $self->variant_file_list,
    );

    unless($result_cmd->execute) {
        die $self->error_message('Failed to create result for model ' . $model->__display_name__ . ' ' . $variant_type . 's using files: ' . join(', ', @files));
    }

    return $result_cmd->manual_result;
}


1;

