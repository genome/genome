package Genome::Model::SomaticValidation::Command::IdentifyDnp;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::IdentifyDnp {
    is => 'Command::V2',
    doc => "Command line interface to create a Genome::Model::Build::SomaticValidation::IdentifyDnpResult",
    has_input => [
        build_id => {
            is => 'Number',
            doc => 'ID of a Somatic Validation Build.',
        },
    ],
    has_output => [
        dnp_result_id => {
            is => 'Text',
            doc => 'Resulting ID of a SoftwareResult for Identify DNP.',
            is_optional => 1,
        },
    ],
    has => [
        dnp_result => {
            is => 'Genome::Model::Build::SomaticValidation::IdentifyDnpResult',
            id_by => 'dnp_result_id',
            doc => "SoftwareResult containing the output of 'gmt somatic identify-dnp'.",
            is_optional => 1,
        },
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
            id_by => 'build_id',
            doc => "Somatic Validation build that the DNP result will belong to and which procides the reads and proportion parameters.",
        },
    ],
    doc => 'identify DNPs',
};

sub sub_command_category { 'pipeline steps' }

sub help_detail {
    "This is a command line wrapper to create a Genome::Model::Build::SomaticValidation::IdentifyDnpResult. It is used by Genome::Model::Build::SomaticValidation's workflow."
}


sub params_for_result {
    my $self = shift;

    my $proportion = $self->build->processing_profile->identify_dnp_proportion;
    unless (defined $proportion) {
        die $self->error_message("'identify_dnp_proportion' not specified on processing profile.");
    }

    my $build = $self->build;
    my @result_users = Genome::SoftwareResult::User->get(user => $build);
    my @results = map($_->software_result, @result_users);
    my @merged_results = grep { $_->isa('Genome::InstrumentData::AlignedBamResult') } @results;
    my ($tumor_alignment_result) = grep( ($_->instrument_data)[0]->sample eq $build->tumor_sample, @merged_results);

    my $tumor_aligment_result_id = $tumor_alignment_result->id;
    my @software_result_users = Genome::SoftwareResult::User->get(user => $self->build);
    my @software_results = map { $_->software_result } @software_result_users;

    my $snv_result_user = Genome::SoftwareResult::User->get(label => 'snv_result', user => $build);
    my $snv_result = $snv_result_user->software_result;

    my $dv2_result_id = $snv_result->id;

    my $test_name = $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME};

    my $params = {
        proportion => $proportion,
        tumor_aligment_result_id => $tumor_aligment_result_id,
        dv2_result_id => $dv2_result_id,
        test_name => $test_name,
    };

    return $params;
}


sub _link_to_result {
    my $self = shift;

    # add user
    my $build = $self->build;
    my $result = $self->dnp_result;
    $result->add_user(user_id => $build->id, user_class_name => $build->class, label => 'uses');

    # symlink result in
    my $identify_dnp_symlink = join('/', $build->data_directory, 'identify-dnp');
    Genome::Sys->create_symlink_and_log_change($self, $result->output_dir, $identify_dnp_symlink);

    return 1;
}


sub shortcut {
    my $self = shift;

    #try to get using the lock in order to wait here in shortcut if another process is creating this alignment result
    my $params = $self->params_for_result;
    my $result = Genome::Model::Build::SomaticValidation::IdentifyDnpResult->get_with_lock(%$params);
    unless($result) {
        $self->debug_message('No existing result found.');
        return;
    }

    $self->debug_message('Using existing result ' . $result->__display_name__);
    $self->dnp_result($result);
    $self->_link_to_result;

    return 1;
}


sub execute {
    my $self = shift;

    unless ($self->build->snv_detection_strategy) {
        $self->warning_message("Snv detection strategy not defined. Skipping IdentifyDnp");
        return 1;
    }

    my $params = $self->params_for_result;
    my $result = Genome::Model::Build::SomaticValidation::IdentifyDnpResult->get_or_create(%$params);

    my $rv = ($result ? 1 : 0);
    unless ($result) {
        $self->error_message('Failed to generate result.');
    }

    $self->dnp_result($result);
    $self->_link_to_result;

    return $rv;
}

1;
