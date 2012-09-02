package Genome::Model::Build::Command::LogTail;

use strict;
use warnings;
use Genome;

class Genome::Model::Build::Command::LogTail {
    is => 'Command::V2',
    has_input => [
        builds          => { is => 'Genome::Model::Build', is_many => 1, doc => "builds with logs to view", shell_args_position => 1 },
    ],
    has_param => [
        n               => { is => 'Number', is_optional => 1, default_value => 10, doc => "the number of lines of tail to show" },
        match_rows      => { is => 'Text', is_optional => 1, default_value => '', doc => 'only show log entries which match this regex' },
        statuses_match  => { is => 'Text', is_optional => 1, default_value => 'crashed|running', doc => 'only show jobs with a status matchin this regex' },
    ],
    doc => 'tail the execution logs of the steps in ' 
};

sub help_detail {
    return shift->help_brief();
}

# TODO: this is pretty much a hack but could be useful if cleaned-up -ss

sub execute {
    my $self = shift;

    my @builds = $self->builds;
    my $n = $self->n;
    my $match_rows = $self->match_rows;
    my $statuses_match = $self->statuses_match;

    $self->status_message("showing rows matching pattern: $match_rows");
    $self->status_message("for steps matching status pattern: $statuses_match");
    $self->status_message("");


    for my $build (@builds) {
        my $prefix;
        if (@builds > 1) {
            $prefix = "BUILD " . $build->__display_name__ . ":\t";
        }
        else {
            $prefix = '';
        }
        
        my $dir = $build->data_directory . '/logs';
        
        # cheap hack
        my @rows = `genome model build view $build->{id}`;
        for my $step (@rows) {
            my ($event_id) = ($step =~ /(\d+)\s*$/);
            next unless $event_id;
            
            print $prefix, $step;
            
            my $edir = $dir . '/' . $event_id . '.err';
            if (-e $edir) {
                print '    ' . $edir . "\n";
                unless ($step =~ /$statuses_match/) {
                    next;
                }
                if ($match_rows ne '') {
                    system "cat $edir | prepend '        ' | egrep '$match_rows' | tail -n $n";
                }
                else {
                    system "tail -n $n $edir | prepend '        '";
                }
            }
        }
        print "\n";
    }
    print "\n";

    return 1;
}

1;

