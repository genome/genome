package Genome::Model::Build::Command::Pids;
use strict;
use warnings;

class Genome::Model::Build::Command::Pids {
    is => 'Command::V2',
    has => [
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
            doc => 'Find pids for these builds',
            shell_args_position => 1,
        },
    ],
    has_optional => [
        strace => {
            is => 'Boolean',
            default => 0,
            doc => 'strace pids - must run as the user of the build. ctrl+c to end one strace and start the next.',
        },
    ],
    doc => 'Print the blades and pids for given builds. Manually run ltrace on the pids for more information.',
};

sub help_brief { 'Print the blades and pids for given builds. Manually run ltrace on the pids for more information.' }
#Automating ltrace seemed to cause some builds to become unstable

sub help_detail{ help_brief() . "\nExample: genome model build pids 10407\nssh BLADE 'ltrace -p PID'" }

sub execute {
    my $self = shift;
    my %seen_host_pids;
    my $dont_trace = "lsf|workflow|grep";

    print "build_id\tjob\thost\n\tpids\n" unless $self->strace;
    for my $build($self->builds){
        my $build_id = $build->id;
        my $build = Genome::Model::Build->get($build_id);
        for my $job ($build->child_lsf_jobs) {#cycle each job id associated with the build
            #There could be multiple execution blades here, but just grab the first one
            #from what I've noticed, multiple blades listed in the execution column are identical anyways
            if(`bjobs -W $job` =~ /edu (?:\d+\*)?(blade.*?)\.gsc.*?(\d+(?:,\d+)+)/){#This regex is environment specific to the Genome Institute
                my $exec_host = $1;
                my @pids = split ',',$2;
                next if $seen_host_pids{$exec_host.$2}++;

                my %pid_desc;
                for my $pid (@pids){
                    my @ps_output = split "\n", `ssh $exec_host 'ps -fp $pid'`;

                    #The end of second line of output from ps is the shell command that executed the process,
                    #which we will use as description. The first line is column headers
                    if($ps_output[1] =~ /(?:\S+\s+){7}(.+)/){
                        $pid_desc{$pid} = $1;
                    }else{
                        $pid_desc{$pid} = '-';
                    }
                }
                my @trace_pids = grep{ $pid_desc{$_} !~ /$dont_trace/i } keys %pid_desc;

                if($self->strace and @trace_pids){
                    system("ssh $exec_host 'strace -p " . join(' -p ', @trace_pids) . "'");
                    print "\n";

                    #Clean up any lingering trace processes that failed to be killed
                    PIDGROUP: for my $pid (@trace_pids){
                        my $grep_output = `ssh $exec_host 'ps aux | grep $pid' | grep strace`;
                        for(split "\n", $grep_output){
                            if(/\S+\s+(\d+).*strace/){
                                `ssh $exec_host 'kill $1'`;
                                last PIDGROUP;
                            }
                        }
                    }
                }elsif($self->strace){
                    print "$build_id\t$job\t$exec_host\n\tno pids\n";
                }else{
                    print "$build_id\t$job\t$exec_host\n";
                    print "\t$_\t$pid_desc{$_}\n" for(@trace_pids);
                }
            } elsif(not $self->strace) {
                print "$build_id\t$job\tno_host\n";
            }
        }
    }

    1;
}
1;
