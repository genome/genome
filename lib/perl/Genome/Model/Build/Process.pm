package Genome::Model::Build::Process;

use strict;
use warnings FATAL => 'all';
use Genome;
use Try::Tiny qw(try catch);

class Genome::Model::Build::Process {
    is => [
        "Genome::Process",
    ],
    has_input => [
        build => {
            is => 'Genome::Model::Build',
        },
    ],

};

sub workflow_name {
    my $self = shift;

    return join(" - ", $self->build->workflow_name, $self->id);
}

sub lsf_project_name {
    my $self = shift;
    return "build" . $self->build->id;
}

sub log_directory {
    my $self = shift;
    return $self->build->log_directory;
}

sub update_status {
    my $self = shift;
    my ($new_status) = @_;

    my $accessor = sprintf("_before_%s", lc($new_status));
    if ($self->can($accessor)) {
        $self->$accessor;
    }

    return $self->SUPER::update_status(@_);
}

sub _before_running {
    my $self = shift;

    try {
        $self->build->initialize;
    } catch {
        $self->_post_build_failure($_);
        die "Couldn't initialize build";
    };

    if ($ENV{WF_USE_FLOW}) {
        $self->build->add_note(
            header_text => 'run using flow',
            body_text => '',
        );
    }
}

sub _before_succeeded {
    my $self = shift;

    unless ( $self->build->success ) {
        my $msg = sprintf(
            'Failed to set build to success: %s',
            $self->build->error_text || 'no error given',
        );
        return $self->_post_build_failure($msg);
    }

    UR::Context->commit() if Genome::Config::get('workflow_builder_backend') ne 'inline';

    require UR::Object::View::Default::Xsl;

    my $cachetrigger = Genome::Config::get('sys_services_web_view_url');
    $cachetrigger =~ s/view$/cachetrigger/;

    my $url = $cachetrigger . '/' . UR::Object::View::Default::Xsl::type_to_url     (ref($self->build)) . '/status.html?id=' . $self->build->id;

    system("curl -k $url >/dev/null 2>/dev/null &");
}

sub _before_crashed {
    my $self = shift;

    my $error = try {
        my $determine_error_command = Genome::Model::Build::Command::DetermineError->create(
            build => $self->build,
            assume_build_status => 'Failed',
            display_results => 0,
            color => 0,
        );
        $determine_error_command->execute();
        return $determine_error_command->formatted_results;
    } catch {
        return "Error could not be determined: $_";
    };
    $self->_post_build_failure($error);
}

sub _post_build_failure {
    my ($self, $msg) = @_;

    $self->error_message($msg);

    my $build_event = $self->build->build_event;
    my $error = Genome::Model::Build::Error->create(
        build_event_id => $build_event->id,
        stage_event_id => $build_event->id,
        stage => 'all stages',
        step_event_id => $build_event->id,
        step => 'main',
        error => $msg,
    );

    unless ( $self->build->fail($error) ) {
        return $self->_failed_build_fail($error);
    }

    return 1;
}

sub _failed_build_fail {
    my ($self, @errors) = @_;

    my $msg = sprintf(
        "Failed to fail build because: %s\nOriginal errors:\n%s",
        ( $self->build->error_text || 'No error given' ),
        join("\n", map { $_->error } @errors),
    );

    $self->error_message($msg);

    return 1;
}

1;
