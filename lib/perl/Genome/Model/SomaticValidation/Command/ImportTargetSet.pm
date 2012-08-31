package Genome::Model::SomaticValidation::Command::ImportTargetSet;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::ImportTargetSet {
    is => 'Command::V2',
    has_input => [
        target_set_file => {
            is => 'Text',
            doc => 'File of the targets from the vendor',
        },
        oid => {
            is => 'Text',
            doc => 'The OID for these targets',
        },
        reference => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'The reference used to make the targets',
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
    doc => 'import the set of targets into the analysis system for validation',
};

sub sub_command_category { 'analyst tools' }

sub help_detail {
    return 'This command is used to upload the target set that is to be used to generate the validation data';
}

sub execute {
    my $self = shift;

    my $name = $self->oid . ' capture chip set';

    my $user = Genome::Sys->username;
    my $sudo_user = Genome::Sys->sudo_username;
    $user .= " ($sudo_user)" if $sudo_user;

   my $description = 'uploaded by ' . $user . ' using ' . $self->command_name;

    my $fl_command = Genome::FeatureList::Command::Create->create(
        name => $name,
        format => 'multi-tracked',
        file_path => $self->target_set_file,
        content_type => 'validation',
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

