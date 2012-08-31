package Genome::Model::Command::Services::Build::Scan;

use strict;
use warnings;

use Guard;
use Workflow;

class Genome::Model::Command::Services::Build::Scan {
    is  => 'Command',
    doc => 'scan all non abandoned, completed or crashed builds for problems',
    has => [
        cron => {
            is => 'Boolean',
            doc =>
              'Reduces output to a subset appropriate for a cron job to email'
        },
        fix => {
            is  => 'Boolean',
            doc => 'Take corrective action in some situations'
        }
    ]
};

sub help_detail {
    return <<EOS
EOS
}

my $printed = 0;

sub write_form {
    if ( !$printed ) {
        print <<'        MARK';
State      Build ID     LSF JOB ID Current Status  Action     Owner    Fix
        MARK
        $printed++;
    }

    write STDOUT;
}

sub execute {
    my $self = shift;

    my (
        $builds_without_job, $builds_with_job,    $job_without_build,
        $events_with_job,    $build_inner_events, $job_info
    ) = $self->build_lists;

    my $state;
    my $build_id;
    my @lsf_id;
    my $event_status;
    my $action;
    my $owner;
    my $fix;

    no warnings;
    format STDOUT =
@<<<<<<<<< @<<<<<<<<<<< @<<<<<<<<< @<<<<<<<<<<<<<< @<<<<<<<<< @<<<<<<< @<<
$state,    $build_id,  shift @lsf_id,   $event_status, $action, $owner, $fix
                        @<<<<<<<<< ~~
                       shift @lsf_id 
.
    use warnings;

    while ( my ( $bid, $lsf_job_ids ) = each %$builds_without_job ) {
        $fix      = ' ';
        $state    = 'no job';
        $build_id = $bid;
        my $event = shift @$lsf_job_ids;
        $event_status = $event->event_status;
        $owner        = $event->user_name;
        my @event_ids_only = keys %$events_with_job;
        my $fixed_status =
          $self->derive_build_status( $build_inner_events->{$build_id},
            \@event_ids_only );
        $action = $self->action_for_derived_status($fixed_status);

        @lsf_id = @$lsf_job_ids;

        if ( $self->fix ) {
            $fix = $self->correct_status( $event->build, $action ) ? 'y' : 'n';
        }
        write_form();
    }

    my $old_builds_with_job  = {};
    my $pend_builds_with_job = {};
    while ( my ( $bid, $lsf_job_ids ) = each %$builds_with_job ) {
        $fix      = ' ';
        $state    = 'run';
        $build_id = $bid;
        my $event = shift @$lsf_job_ids;
        $event_status = $event->event_status;
        $owner        = $event->user_name;
        $action       = 'none';
        my %joined = map { $_ => 1 } @{ $lsf_job_ids->[0] },
          @{ $lsf_job_ids->[1] };

        @lsf_id = sort {
            exists $job_info->{$b} <=> exists $job_info->{$a}
              || $job_info->{$a}->started_on <=> $job_info->{$b}->started_on
        } keys %joined;

        my $t    = $job_info->{ $lsf_id[0] }->started_on;
        my $ftd  = 60 * 60 * 24 * 15;
        my $ftda = time - $ftd;
        if ( !defined $t ) {
            $pend_builds_with_job->{$bid} = [ $event, @lsf_id ];
        } elsif ( $t < $ftda ) {
            $old_builds_with_job->{$bid} = [ $event, @lsf_id ];
        } else {
            write_form() unless $self->cron;
        }
    }

    while ( my ( $bid, $lsf_job_ids ) = each %$pend_builds_with_job ) {
        $fix      = ' ';
        $state    = 'pend';
        $build_id = $bid;
        my $event = shift @$lsf_job_ids;
        $event_status = $event->event_status;
        $owner        = $event->user_name;
        $action       = 'none';
        @lsf_id       = @$lsf_job_ids;

        write_form() unless $self->cron;
    }

    while ( my ( $bid, $lsf_job_ids ) = each %$old_builds_with_job ) {
        $fix      = ' ';
        $state    = 'old';
        $build_id = $bid;
        my $event = shift @$lsf_job_ids;
        $event_status = $event->event_status;
        $owner        = $event->user_name;
        $action       = 'kill';
        @lsf_id       = @$lsf_job_ids;

        write_form();
    }

    while ( my ( $bid, $lsf_job_ids ) = each %$job_without_build ) {
        $fix      = ' ';
        $state    = 'no build';
        $build_id = $bid;
        my $event = shift @$lsf_job_ids;
        if (defined $event) {
            $event_status = $event->event_status;
            $owner        = $event->user_name;
        } else {
            $event_status = '(undef)';
            $owner = '(undef)';
        }
        my @event_ids_only = keys %$events_with_job;
        my $fixed_status =
          $self->derive_build_status( $build_inner_events->{$build_id},
            \@event_ids_only );
        $action = 'kill';
        my $daction = $self->action_for_derived_status($fixed_status);

        if ( $daction ne $self->action_for_derived_status($event_status) ) {
            $action .= ' ' . $daction;
        }
        @lsf_id = @$lsf_job_ids;

        if ( $self->fix ) {
            if ($action eq 'kill fail') {
                my $ok = scalar @lsf_id;
                foreach my $job_id (@lsf_id) {
                    my $rv = system("bkill $job_id");
                    if (!$rv) {
                        $ok--;
                    }
                }
                $fix = $ok == 0 ? 'y' : 'n';
            }
        }

        write_form();
    }

    return 1;
}

sub correct_status {
    my ( $self, $build, $action ) = @_;

    #    if (!defined $action) {
    #        return 0;
    if ( $action eq 'success' ) {
        $build->success;
        return 1;
    } elsif ( $action eq 'fail' ) {
        $build->fail;
        return 1;
    } elsif ( $action eq 'abandon' ) {
        $build->build_event->event_status('Abandoned');
        $build->build_event->date_completed( UR::Time->now );
        return 1;
    }
    return 0;
}

sub action_for_derived_status {
    my $self = shift;
    my $derived = shift || return 'none';

    if (   $derived eq 'Failed'
        || $derived eq 'Crashed'
        || $derived eq 'Running'
        || $derived eq 'Scheduled'
        || $derived eq 'unknown' )
    {
        return 'fail';
    } elsif ( $derived eq 'Succeeded' ) {
        return 'success';
    } elsif ( $derived eq 'Abandoned' ) {
        return 'abandon';
    } else {
        return 'none';
    }
}

sub build_lists {
    my $self = shift;

    my $builds_db = {};
    {
        my $iter = Genome::Model::Event->create_iterator(
            event_type   => 'genome model build',
            event_status => [ 'Scheduled', 'Running' ]
        );
        while ( my $event = $iter->next ) {
            $builds_db->{ $event->build_id } = [ $event->lsf_job_id, $event ];
        }
    }

    my $events_lsf = {};
    my $builds_lsf = {};
    my $jobs_lsf   = {};
    {
        my $ji = Job::Iterator->new;
        while ( my $job = $ji->next ) {
            next if ( $job->{Status} eq 'PSUSP' );
            if ( $job->{Command} =~ /build run.+?--build-id (\d+)/ ) {
                $builds_lsf->{$1} ||= [];
                push @{ $builds_lsf->{$1} }, $job->{Job};
            }
            if ( exists $job->{'Job Name'} && $job->{'Job Name'} =~ /(\d+)$/ ) {
                $events_lsf->{$1} ||= [];
                push @{ $events_lsf->{$1} }, $job->{Job};
            }

            $jobs_lsf->{ $job->{Job} } = $job;
        }

        #        open my $bjobs, "bjobs -u all -w |";
        #        while ( my $line = <$bjobs> ) {
        #            if ( $line =~ /^(\d+).+?build run.+?--build-id (\d+)/ ) {
        #                $builds_lsf->{$2} ||= [];
        #                push @{ $builds_lsf->{$2} }, $1;
        #            }
        #            if ( $line =~
        #/^(\d+).+?(\d+) (?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)/
        #              )
        #            {

        # this is a terrible match, but its the best i can do right now
        #                $events_lsf->{$2} ||= [];
        #                push @{ $events_lsf->{$2} }, $1;
        #            }
        #        }
        #        close $bjobs;
    }

    my $builds_both = { map { $_ => 1 } keys %$builds_lsf, keys %$builds_db };

    my $inner_events_db = {};
    {
        my $iter = Genome::Model::Event->create_iterator(
            build_id   => [ keys %$builds_both ],
            event_type => { operator => 'ne', value => 'genome model build' }
        );
        while ( my $event = $iter->next ) {
            $inner_events_db->{ $event->build_id } ||= {};
            $inner_events_db->{ $event->build_id }->{ $event->id } = $event;
        }
    }

    foreach my $bid ( keys %$builds_both ) {
        my $event =
          exists $builds_db->{$bid}
          ? $builds_db->{$bid}->[1]
          : Genome::Model::Event->get(
            event_type => 'genome model build',
            build_id   => $bid
          );
        if ( exists $builds_lsf->{$bid} && exists $builds_db->{$bid} ) {
            my $from_lsf = delete $builds_lsf->{$bid};
            my $from_db  = delete $builds_db->{$bid};
            pop @$from_db;

            $builds_both->{$bid} = [ $event, $from_lsf, $from_db ];
        } else {

            if ( exists $builds_lsf->{$bid} ) {
                unshift @{ $builds_lsf->{$bid} }, $event;
            } elsif ( exists $builds_db->{$bid} ) {
                unshift @{ $builds_db->{$bid} }, $event;
                pop @{ $builds_db->{$bid} };
            }
            delete $builds_both->{$bid};
        }
    }

    return $builds_db, $builds_both, $builds_lsf, $events_lsf, $inner_events_db,
      $jobs_lsf;
}

sub derive_build_status {
    my $self              = shift;
    my $events_in_build   = shift;
    my $running_event_ids = shift;

    my %status_prio = (
        Failed    => 90,
        Crashed   => 80,
        Running   => 70,
        Scheduled => 30,
        Succeeded => 20,
        Abandoned => 10,
        unknown   => -100,
    );

    # Using the above status priority list we will derive a new
    # build_event status.  At first glance it seems like Abandoned
    # is too low, but we only want the whole build to show that
    # status if there is nothing else higher.  In the past it was
    # common to Abandon a single lane and restart the build

    my $new_status = 'unknown';
    while ( my ( $eid, $event ) = each %$events_in_build ) {
        my $event_status = $event->event_status;
        if ( $event_status eq 'Running' && defined $running_event_ids ) {
            ## see if its really running
            my $found = 0;
            foreach (@$running_event_ids) {
                if ( $_ == $event->id ) {
                    $found = 1;
                    last;
                }
            }

            $event_status = 'Crashed' if !$found;
        }

        if ( exists $status_prio{$event_status}
            && $status_prio{$event_status} > $status_prio{$new_status} )
        {
            $new_status = $event_status;
        }
    }
    return $new_status;
}

###
### The following will be moved into their own files (and possibly renamed) if
### they prove to be more stable than the other parsing code used in GME and
### Workflow::Server::Hub

package Job::Iterator;

sub new {
    my $class = shift;
    my $opts = shift || '-u all';
    open my $f, 'bjobs -l ' . $opts . ' |';
    readline $f;    ## read off the blank line then read the first real line
    return bless [ $f, scalar( readline($f) ) ], $class;
}

sub next {
    my $self = shift;

    while (1) {
        my $line = readline( $self->[0] );
        push @$self, $line;

        if ( defined $line && $line =~ /^Job \</ ) {
            last;
        } elsif ( !defined $line ) {
            last;
        }
    }
    my @jl = splice( @$self, 1, -1 );
    if ( defined $jl[0] ) {
        return Job->new(@jl);
    } else {
        return;
    }
}

package Job;

sub new {
    my ( $class, @lines ) = @_;
    my $spool = join( '', @lines );

    # this regex nukes the indentation and line feed
    $spool =~ s/\s{22}//gm;

    my @eventlines = split( /\n/, $spool );
    my %jobinfo = ();

    my $jobinfoline = shift @eventlines;
    if ( defined $jobinfoline ) {

        # sometimes the prior regex nukes the white space between Key <Value>
        $jobinfoline =~ s/(?<!\s{1})</ </g;
        $jobinfoline =~ s/>,(?!\s{1})/>, /g;

        # parse out a line such as
        # Key <Value>, Key <Value>, Key <Value>
        while ( $jobinfoline =~
            /(?:^|(?<=,\s{1}))(.+?)(?:\s+<(.*?)>)?(?=(?:$|;|,))/g )
        {
            $jobinfo{$1} = $2;
        }
    }

    $jobinfo{__events} = [];
    foreach my $el (@eventlines) {
        if ( $el =~
/^(Sun|Mon|Tue|Wed|Thu|Fri|Sat) (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s{1,2}(\d{1,2})\s{1,2}(\d{1,2}):(\d{2}):(\d{2}):/
          )
        {

            $el =~ s/(?<!\s{1})</ </g;
            $el =~ s/>,(?!\s{1})/>, /g;

            my $time = substr( $el, 0, 21, '' );
            substr( $time, -2, 2, '' );

            # see if we really got the time string
            if ( $time !~ /\w{3} \w{3}\s+\d{1,2}\s+\d{1,2}:\d{2}:\d{2}/ ) {

                # there's stuff we dont care about at the bottom, just skip it
                next;
            }

            my $desc = {};
            while ( $el =~
/(?:^|(?<=,\s{1}))(.+?)(?:\s+<(.*?)>(?:\s+(.*?))?)?(?=(?:$|;|,))/g
              )
            {
                if ( defined $3 ) {
                    $desc->{$1} = [ $2, $3 ];
                } else {
                    $desc->{$1} = $2;
                }
            }
            push @{ $jobinfo{__events} }, Job::Event->new( $time, $desc );

        }
    }

    return bless \%jobinfo, $class;
}

sub started_on {
    my $self = shift;

    for ( reverse @{ $self->{__events} } ) {
        if ( exists $_->{'Started on'} ) {
            return $_->{__time};
        }
    }
    return;
}

sub events {
    @{ $_[0]->{__events} };
}

package Job::Event;

use Date::Parse;

sub new {
    my ( $class, $time, $desc ) = @_;

    $desc->{__time} = str2time($time);
    return bless $desc, $class;
}

1;
