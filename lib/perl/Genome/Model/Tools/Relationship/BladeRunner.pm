package Genome::Model::Tools::Relationship::BladeRunner;

use strict;
use warnings;
use Data::Dumper;
use Genome;           
use Genome::Info::IUB;
use POSIX;
use Cwd;
use File::Basename;
use File::Path;
use LSF::Job;

class Genome::Model::Tools::Relationship::BladeRunner {
    is => 'Command',
    has => [
    jobs => {
     },
    running_jobs=> {
        is_optional=>1,
    },
    finished_jobs => {
        is_optional=>1,
    },
    queue=> {
        default=>"long",
    },
    ],
};

sub help_brief {
}

sub help_detail {
}


sub execute {
    my $self = shift;
    my ($jobs_ready_to_start, $dependent_jobs) = $self->find_jobs_to_start($self->jobs);
    
    while(scalar(@$jobs_ready_to_start) + scalar(@$dependent_jobs) > 0) { ##main event loop
        for my $job (@{$jobs_ready_to_start}) {
            $self->error_message("Starting a job... " . $job->command_line);
            my $job = $job->execute;
            if(!defined($job)) { # unable to submit for some reason
                $self->error_message("Job Unable to start:\n" . Dumper($job));
                return 0;
            }
            else {
                my @jobs;
                if($self->running_jobs) {
                    @jobs = @{$self->running_jobs};
                }
                push @jobs, $job;
                $self->running_jobs(\@jobs);
            }
        }
        while(1) {
            $self->error_message("Entering sleep/event loop...");
            sleep(10);
            my $new_finished_jobs = $self->find_finished_jobs($self->finished_jobs, $self->running_jobs);
            if($new_finished_jobs) {
                $self->error_message("Found successful finished jobs. Leaving sleep loop to look for more jobs to launch");
                ($jobs_ready_to_start, $dependent_jobs) = $self->find_jobs_to_start($dependent_jobs);
                last;
            }
        }
    }
    return 1;
}

sub find_finished_jobs {
    my ($self, $finished_jobs, $running_jobs) = @_;
    ##confirm the jobs are succeeded or else we have to die here
    $DB::single=1; 
    my @ids_of_interest = map {$_->job_id} @$running_jobs;
    my $bjobs_cmd = "bjobs " . join(" ", @ids_of_interest);
    my @bjobs_return = `$bjobs_cmd`;
    my ($new_finished_jobs, $new_running_jobs) = $self->parse_bjobs(\@bjobs_return, $running_jobs);
    if(scalar(@$new_finished_jobs) > 0) {
            $self->error_message("Found finished jobs to check for success..." . scalar(@$new_finished_jobs));
        my $all_successful=$self->confirm_success($new_finished_jobs);
        if(!$all_successful) {
            $self->error_message("Some Jobs failed. Bailing out\n");
            die;
        }
        $self->finished_jobs($new_finished_jobs);
        $self->running_jobs($new_running_jobs);
    }
    return scalar(@$new_finished_jobs) ;
}
    
sub confirm_success {
    my ($self, $new_finished_jobs) = @_;
    for my $job (@$new_finished_jobs) {
        if(!$job->bjob_l_is_job_successful) {
            $self->error_message("Job has failed\n" . Dumper($job));
            return 0;
        }
    }
    return 1;
}


sub parse_bjobs {
    my ($self, $bjobs_return, $running_jobs) = @_;
    my %finished_job_ids;
    $DB::single=1;
    for my $line (@$bjobs_return) {
        my ($id, $user_name, $status, @rest) = split /\s+/, $line;
        if($status eq 'DONE' || $status eq 'EXIT') {
            $finished_job_ids{$id}=1;
        }
    }
    my (@new_finished_jobs, @still_running_jobs);
    for my $job (@$running_jobs) {
        if(exists($finished_job_ids{$job->job_id})) {
            push @new_finished_jobs, $job;
        }
        else {
            push @still_running_jobs, $job;
        }
    }
    return (\@new_finished_jobs, \@still_running_jobs);
}


sub find_jobs_to_start {
    my ($self, $jobs)= @_;
    my (@ready_to_start_jobs, @dependent_jobs);
    for my $job (@{$jobs}) {
        if(!defined($job->parents)) {
            push @ready_to_start_jobs, $job;
        }
        elsif($self->parents_are_done($job)) { #job has parents
            push @ready_to_start_jobs, $job;
        }
        else {
            push @dependent_jobs, $job;
        }
    }
    return (\@ready_to_start_jobs, \@dependent_jobs);
}

sub parents_are_done { 
    my ($self, $job) = @_;
    my $finished_jobs = $self->finished_jobs;
    my @parent_jobs = @{$job->parents};
    my %jobs_that_need_to_finish;
    for my $parent_job (@parent_jobs) {
        $jobs_that_need_to_finish{$parent_job->command_line}=1;
    }
    for my $finished_job (@$finished_jobs) {
        if(exists($jobs_that_need_to_finish{$finished_job->command_line})) {
            delete $jobs_that_need_to_finish{$finished_job->command_line};
        }
    }
    if(scalar(keys %jobs_that_need_to_finish)==0) {    
        return 1;
    }
    return 0;
}

