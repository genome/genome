package Genome::Model::Command::Services::Build::Scan;

use strict;
use warnings;
use Genome;
use IPC::System::Simple qw(capture);

class Genome::Model::Command::Services::Build::Scan {
    is => 'Command::V2',
    has => [
        fix => {
            is => 'Boolean',
            doc => 'Take corrective action in some situations'
        },
    ],
};

sub execute {
    my $self = shift;
    my $running_builds = $self->_get_running_builds;
    my $failed_build_count = 0;
    my $failed_to_fail_count = 0;
    my $should_fail_count = 0;
    print join("\t", qw(Build status run_by date_scheduled action))."\n";
    for my $build (@$running_builds) {
        my $should_fail = 0;
        eval {
            $should_fail = $self->_build_is_bad($build);
            if ($should_fail){
                $should_fail_count++;
                print join("\t", $build->id, $build->status, $build->run_by, $build->date_scheduled, "fail")."\n";
                if ($self->fix) {
                    eval {$build->fail};
                    if ($build->status eq "Failed") {
                        $failed_build_count++;
                    }
                    else {
                        $failed_to_fail_count++;
                        die sprintf("Problem failing build %s\n", $build->id);
                    }
                }
            }
        };
        if ($@) {
            print sprintf("Problem with build %s: %s\n", $build->id, $@);
        }
    }
    print sprintf("Should fail %d builds\nSuccessfully failed %d builds\nFailed to fail %d builds\n",
            $should_fail_count, $failed_build_count, $failed_to_fail_count);
    return 1;
}

sub _build_is_bad {
    my $self = shift;
    my $build = shift;
    my $has_jobs = $self->_has_lsf_jobs($build);
    if ($has_jobs) {
        return 0;
    }
    #Build could have succeeded in between getting the list of running builds and
    #checking the jobs
    UR::Context->current->reload($build->the_master_event);
    if ($build->status eq "Scheduled" or $build->status eq "Running") {
        return 1;
    }
    else {
        return 0;
    }
}

sub _get_running_builds {
    my $self = shift;
    my @builds = Genome::Model::Build->get(status => [qw(Running Scheduled)]);
    return \@builds;
}

sub _has_lsf_jobs {
    my $self = shift;
    my $build = shift;
    my @output = eval {capture('bjobs', '-uall', '-Pbuild'.$build->id)};
    if ($@) {
        die $self->error_message("Problem running bjobs output");
    }
    if (@output) {
        return 1;
    }
    else {
        return 0;
    }
}

1;

