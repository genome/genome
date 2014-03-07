package Genome::Sys::LSF::JobIterator;

use strict;
use warnings;

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
            /^(Sun|Mon|Tue|Wed|Thu|Fri|Sat) (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s{1,2}(\d{1,2})\s{1,   2}(\d{1,2}):(\d{2}):(\d{2}):/
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

