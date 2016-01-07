package Genome::Test::Command::RunModelTest;

use strict;
use warnings;
use Data::Dump qw(pp);

use Genome;

class Genome::Test::Command::RunModelTest {
    is => 'Command::V2',
    roles => 'Genome::Role::CommandWithColor',
    doc => "Run a single model test",
    has => [
        model => {
            is => 'Genome::Model',
            require_user_verify => 0,
            doc => 'Model to test',
            shell_args_position => 1,
        },
        timeout => {
            is => 'Number',
            default => 120,
            doc => 'Number of minutes to wait on the build before failing.',
        },
        polling_interval => {
            is => 'Number',
            default => 30,
            doc => 'Number of seconds to sleep between polling for '.
                'build status.',
        },
        software_revision => {
            is => 'Text',
            is_optional => 1,
            doc => 'Used similiar to GENOME_SOFTWARE_RESULT_TEST_NAME '.
                'but for builds instead of SoftwareResults.',
        },
    ],
};

sub execute {
    my $self = shift;

    # flush output buffer after every write or print
    local $| = 1;

    $self->print_inputs();

    my $build = $self->get_or_create_build();

    $self->wait_for_build_to_complete($build);

    if ($build->status eq 'Succeeded') {
        $self->diff_with_blessed_build($build);
    } else {
        die sprintf("Build (%s) was not successful",
            $build->id);
    }

    return 1;
}

sub print_inputs {
    my $self = shift;

    $self->print_input("model", $self->model->__display_name__);
    for my $name qw(timeout polling_interval software_revision) {
        $self->print_input($name);
    }
    print "\n";
}

sub print_input {
    my ($self, $name, $value) = @_;
    $value = $self->$name unless defined($value);
    printf "%s\n", $self->_color_pair($name, pp($value));
}

sub get_or_create_build {
    my $self = shift;

    my $build;
    if (defined($self->software_revision)) {
        $build = $self->get_existing_build();
    }

    if (defined($build)) {
        printf $self->_color("Found existing suitable build: %s\n\n", 'cyan'),
            pp($build->__display_name__);
        return $build;
    } else {
        printf "No existing suitable build found\n\n";
        return $self->create_build();
    }
}

sub create_build {
    my $self = shift;

    my %params = (
        model_id => $self->model->id,
        software_revision => $self->software_revision,
    );
    printf $self->_color("Creating build with params: %s\n", 'green'),
        pp(\%params);
    my $build = Genome::Model::Build->create(%params);
    die "Could not create new build!" unless ($build);

    print("Starting build\n");
    unless ($build->start()) {
        die "Cound not start new build!";
    }

    print("Saving build\n");
    unless (UR::Context->commit()) {
        die "Could not save new build!";
    }

    return $build;
}

sub get_existing_build {
    my $self = shift;

    my %params = (
        model_name => $self->model->name,
        run_by => Genome::Sys->username,
        software_revision => $self->software_revision,
        'status in' => ['Scheduled', 'Running', 'Succeeded'],
        '-order' => '-date_completed',
    );
    printf "Looking for build with params: %s\n", pp(\%params);

    my @builds = Genome::Model::Build->get(%params);
    if (scalar(@builds)) {
        if (scalar(@builds) > 1) {
            printf "Found %s builds, choosing the most recently completed one\n",
                scalar(@builds);
        }
        return $builds[0];
    } else {
        return;
    }
}

sub wait_for_build_to_complete {
    my $self = shift;
    my $build = shift;

    my $start_time = time();
    my $timeout_seconds = $self->timeout * 60.0;

    my $event = $build->the_master_event;
    unless ($event) {
        die "Could not get the build's master event!";
    }

    printf "Monitoring build (%s) until it completes or timeout "
        . "of %s minutes is reached.", $build->id, $self->timeout;

    while (!grep { $event->event_status eq $_ } ('Succeeded',
            'Failed', 'Crashed')) {
        UR::Context->current->reload($event);
        UR::Context->current->reload($build);
        my $elapsed_time = time() - $start_time;
        if ($elapsed_time > $timeout_seconds) {
            die sprintf("Build (%s) timed out after %s minutes",
                $build->id, $self->timeout);
        }

        print '.';
        sleep($self->polling_interval);
    }
    print "\n";
}

sub diff_with_blessed_build {
    my ($self, $build)= @_;

    my $diff_cmd = Genome::Model::Build::Command::DiffBlessed->create(
        new_build => $build,
    );
    unless ($diff_cmd->execute) {
        die "Diff command failed to execute!";
    }
    die "Diffs found, exiting." if $diff_cmd->has_diffs;
}
