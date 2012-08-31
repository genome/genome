package Genome::Model::SomaticValidation::Command::ImportDesign;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::ImportDesign {
    is => 'Command::V2',
    has => [
        design_file => {
            is => 'Text',
            doc => 'BED file of designs sent to the vendor',
        },
        project_name => {
            is => 'Text',
            doc => 'Name of the project (will be used to name the design in the system)',
        },
        reference => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'The reference used to make the design',
        },
    ],
    has_optional_output => [
        feature_list_id => {
            is => 'Text',
            doc => 'ID of the newly created feature-list',
        },
        feature_list => {
            is => 'Genome::FeatureList',
            id_by => 'feature_list_id',
            doc => 'the newly created feature-list',
        },
    ],
    doc => 'import the design into the analysis system for validation',
};

sub sub_command_category { 'analyst tools' }

sub help_detail {
    return 'This command is used to upload the design sent off to prepare for validation';
}

sub execute {
    my $self = shift;

    my $name = $self->project_name . ' capture design';

    my $user = Genome::Sys->username;
    my $sudo_user = Genome::Sys->sudo_username;
    $user .= " ($sudo_user)" if $sudo_user;

    my $description = 'uploaded by ' . $user . ' using ' . $self->command_name;

    my $fl_command = Genome::FeatureList::Command::Create->create(
        name => $name,
        format => 'true-BED',
        file_path => $self->design_file,
        content_type => 'roi',
        reference => $self->reference,
        description => $description,
    );

    $fl_command->dump_status_messages(1);
    $fl_command->dump_warning_messages(1);
    $fl_command->dump_error_messages(1);

    my $fl;
    unless($fl = $fl_command->execute) {
        die $self->error_message('Failed to import target set file.');
    }

    $self->feature_list($fl);
    return 1;
}

1;

