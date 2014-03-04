package Genome::Model::Command::Services::Build::Scan;

use strict;
use warnings;
use Genome;

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
        unless ($self->_has_lsf_jobs($build)){
            $should_fail_count++;
            print join("\t", $build->id, $build->status, $build->run_by, $build->date_scheduled, "fail")."\n";
            if ($self->fix) {
                eval {
                    $build->fail;
                };
                if ($build->status eq "Failed") {
                    print sprintf("Failed build %s\n", $build->id);
                    $failed_build_count++;
                }
                else {
                    print sprintf("Problem failing build %s\n", $build->id);
                    $failed_to_fail_count++;
                }
            }
        }
    }
    print sprintf("Should fail %d builds\nSuccessfully failed %d builds\nFailed to fail %d builds\n",
            $should_fail_count, $failed_build_count, $failed_to_fail_count);
    return 1;
}

sub _get_running_builds {
    my $self = shift;
    my @builds = Genome::Model::Build->get(status => [qw(Running Scheduled)]);
    return \@builds;
}

sub _has_lsf_jobs {
    my $self = shift;
    my $build = shift;
    my $cmd = 'bjobs -u all -P build'.$build->id. ' 2>&1';
    my $output = `$cmd`;
    if ($output =~ /^JOBID/) {
        return 1;
    }
    elsif ($output =~ /^No job found/) {
        return 0;
    }
    else {
        print "Output was:$output\n";
        die $self->error_message("Problem parsing bjobs output for build ".$build->id);
    }
}

1;

