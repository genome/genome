package Genome::Model::ReferenceSequence::Command::CreateFeatureListInput;

use strict;
use warnings;
use Genome;
use Set::Scalar;
use IPC::Run qw(run);

class Genome::Model::ReferenceSequence::Command::CreateFeatureListInput {
    is => 'Command::V2',
    has => [
        bed_file => {
            is => 'Path',
            doc => 'The (true-BED) bed file path that should be made into a feature list. This should not have a header.',
        },
        feature_list_name => {
            is => 'String',
            doc => 'The name of the feature list to be created',
        },
        build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'The build to which the feature list should be assigned',
        },
        input_name => {
            is => 'String',
            doc => 'The name of the input property on the build to be set',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->validate;

    my $clean_bed = $self->cleanup_bed;

    my $feature_list_command = Genome::FeatureList::Command::Create->create(
        name => $self->feature_list_name,
        format => 'true-BED',
        content_type => 'roi',
        file_path => $clean_bed,
        reference => $self->build,
    );

    my $feature_list = $feature_list_command->execute;

    $self->validate_feature_list($feature_list);

    my $input_name = $self->input_name;
    $self->build->$input_name($feature_list);

    return 1;
}

sub validate {
    my $self = shift;

    my $input_name = $self->input_name;
    unless ($self->build->can($input_name)) {
        die $self->error_message("The build class provided (%s) has no input name (%s)", $self->build->class, $input_name);
    }
    Genome::Sys->validate_file_for_reading($self->bed_file);

    return 1;
}

sub validate_feature_list {
    my ($self, $feature_list) = @_;

    my $feature_list_chromosomes = Set::Scalar->new($feature_list->chromosome_list);
    my $reference_chromosomes = Set::Scalar->new(@{$self->build->chromosome_array_ref});
    $DB::single=1;
    unless ($feature_list_chromosomes <= $reference_chromosomes) {
        die $self->error_message("Some chromosomes exist in the feature list and are not present in the reference: %s", $feature_list_chromosomes - $reference_chromosomes);
    }

    return 1;
}

sub cleanup_bed {
    my $self = shift;

    my $clean_bed_path = Genome::Sys->create_temp_file_path;
    unless (run(['sed', 's/^chr//', $self->bed_file], '>', $clean_bed_path)) {
        die $self->error_message("sed failure: %s", $?);
    }
    return $clean_bed_path;
}


1;

