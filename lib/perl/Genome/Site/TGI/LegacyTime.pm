# A few useful time functions.
 
package UR::Time;

=pod

=head1 NAME

UR::Time - a few useful time-related functions.

=head1 SYNOPSIS

  ##- use UR::Time;

  UR::Time->config(time => '%H:$M:%S');
  $format = UR::Time->config('time');

  $date = UR::Time->today;

  $date_time = UR::Time->now;

  $ds = UR::Time->timediff($t1, $t2);

  $rv = UR::Time->compare_dates($d1, $d2);

  ($sec, $min, $hour, $day, $mon, $year)
      = UR::Time->datetime_to_numbers($datetime);

  $datetime = UR::Time->numbers_to_datetime($sec, $min, $hour,
                                             $day, $mon, $year);

=head1 DESCRIPTION

This package provides several useful functions for reporting and
comparing dates and times in a standard date format.  The format is
set using the C<config> method.  All years are expected to be not
truncated and all times are expected in 24 hour format.

=cut

# set up package
require 5.006_000;
use warnings;
use strict;
our $VERSION = "0.29"; # UR $VERSION;

# set up module
use UR;
use base qw(UR::ModuleBase UR::ModuleConfig);
use Date::Pcalc;
use Date::Parse;
use POSIX;
use Time::Local;

=pod

=head1 METHODS

These methods use Perl builtin time methods and C<Date::Pcalc> (see
L<Date::Pcalc>) to get date and time information and then format it or
operate on it.

=over 4

=item config

  UR::Time->config(datetime => '%Y-%m-%d %H:%M:%S');
  $dt_format = UR::Time->config('datetime');
  %config = UR::Time->config;

This method is used to set and retrieve configuration information.
See L<UR::ModuleConfig> for details of method behavior.  Possible
configuration options are:

=over 6

=item datetime

The format for date/time represenations in scalar format.  The default
is shown above.

=item date

The format for date represenations in scalar format.  The default is
C<%Y-%m-%d>.

=item time

The format for time represenations in scalar format.  The default is
C<%H:%M:%S>.

=back

All formats are interpreted as format arguments to strftime (see
L<strftime>).

=cut

# set configuration variables
__PACKAGE__->config
(
    datetime => '%Y-%m-%d %H:%M:%S',
    date => '%Y-%m-%d',
    time => '%H:%M:%S',
);

=pod

=item today

  $date = UR::Time->date_now;

This method returns today's date in the default format.
(Aliassed as "today".)

=cut

# returns today's date in default Oracle date format
*today = \&date_now;
sub date_now
{
    my $class = shift;

    # get date format
    my $format = $class->config('date');
    if ($format)
    {
        $class->debug_message("got date format %s", $format);
    }
    else
    {
        $class->error_message("no date format specified");
        return;
    }

    # return formatted current time
    return POSIX::strftime($format, localtime);
}

=pod

=item time_now

  $date = UR::Time->time_now;

Returns the time "now", much like the now method w/o the date.

=cut

# returns today's date in default Oracle date format
sub time_now
{
    my $class = shift;

    # get date format
    my $format = $class->config('time');
    if ($format)
    {
        $class->debug_message("got date format %s", $format);
    }
    else
    {
        $class->error_message("no date format specified");
        return;
    }

    # return formatted current time
    return POSIX::strftime($format, localtime);
}


=pod

=item now

  $date = UR::Time->now;

This method returns the current date/time in the default format.  This
method will attempt to get this time from the the database, if a
connection is available.

=cut

# returns current date and time in default Oracle date format
my $use_local_time = 0;
my @local_time_offset;
sub now
{
    my $class = shift;

    # do not consult database if we already have
    if ($use_local_time) 
    {   
	# calculate the time from the calculated offset
	my ($y, $mon, $d, $h, $min, $s)
            = Date::Pcalc::Add_Delta_YMDHMS(Date::Pcalc::Today_and_Now,
                                           @local_time_offset);
        # format the date accordingly
        return $class->numbers_to_datetime($s, $min, $h, $d, $mon, $y);
    }
    # else just get the time from the machine
    return $class->now_local;
}

=pod

=item now_local

  $date = UR::Time->now_local;

This method returns the current machine date/time in the default
format.  It never consults the database and will match current local
file timestamps.

=cut

# returns current date and time in default Oracle date format
sub now_local
{
    my $class = shift;
    return $class->from_integer(time);
}

=pod

=item from_integer

  $date = UR::Time->from_integer($time);

Takes a time in long-integer format (as returned by stat() or time()),
and returns the date/time in the default format.

=cut

# returns epoch seconds in default datetime format
sub from_integer
{
    my $class = shift;
    my ($time) = @_;

    # get date format
    my $format = $class->config('datetime');
    if ($format)
    {
        $class->debug_message("got date format %s", $format);
    }
    else
    {
        $class->error_message("no datetime format specified");
        return;
    }

    # return formatted current time
    return POSIX::strftime($format, localtime($time));
}

=pod

=item timediff

  $ds = UR::Time->timediff($t1, $t2);

This method computes the difference in seconds between two times
(C<$t1> and C<$t2>) in 24 hour format (C<HH:MM:SS>).  The difference
in seconds is returned.  If the times are not in the correct format,
C<undef> is returned.  If C<$t1> is greater than C<$t2>, a negative
value is returned.

=cut

# return difference in seconds of two times
sub timediff
{
    my $class = shift;
    my ($t1, $t2) = @_;

    foreach my $t ($t1, $t2)
    {
	if ($t =~ /(\d\d):(\d\d):(\d\d)/)
	{
	    $t = $3 + ($2 * 60) + ($1 * 60 * 60);
	}
	else
	{
            $class->error_message("failed to parse times:$t1,$t2");
	    return;
	}
    }

    return $t2 - $t1;
}

=pod

=item compare_dates

  $rv = UR::Time->compare_dates($d1, $d2);

This method compares two date/times and returns -1 if C<$d1> is less
than C<$d2>, 1 if C<$d1> is greater than C<$d2>, and 0 if the dates
are equivalent.

=cut

# returns equivalent of first date <=> second date
# an undef date is considered to be earlier than a defined date
sub compare_dates
{
    my $class = shift;
    my ($date1, $date2) = @_;

    no warnings;
    $class->debug_message("comparing >$date1< with >$date2<");
    use warnings;
    
    return 0 if !$date1 && !$date2;
    # undefined dates are less than defined ones
    return -1 if !$date1;
    return 1 if !$date2;

    # split up the date into parts
    my ($s1, $n1, $h1, $d1, $m1, $y1) = UR::Time->datetime_to_numbers($date1);
    my ($s2, $n2, $h2, $d2, $m2, $y2) = UR::Time->datetime_to_numbers($date2);
    return undef unless $s1 && $s2; # error

    $class->debug_message("$y1,$m1,$d1,$h1,$n1,$s1 to $y2,$m2,$d2,$h2,$n2,$s2");

    no warnings qw(uninitialized);
    return $y1 <=> $y2 || $m1 <=> $m2 || $d1 <=> $d2
	|| $h1 <=> $h2 || $n1 <=> $n2 || $s1 <=> $s2;
}

=pod

=item datetime_to_numbers

  ($sec, $min, $hour, $day, $mon, $year)
      = UR::Time->datetime_to_numbers($datetime);

This method converts a date/time in the default format into a list of
the numeric represenation of the year, month, day, hour, minute, and
second.  If any part of the date or time cannot be determined from the
string, it's value will be undef.  All values start at 1, e.g.,
January is month 1, not zero.

=cut

# split up a date in the default Oracle date format into parts
# return numeric year, month, day, hour, minute, and second
sub datetime_to_numbers
{
    my $class = shift;
    my ($date) = @_;

    # try to parse
    my @time = strptime($date);
    if (@time)
    {
        $class->debug_message("parse date $date");
        # fix month
        ++$time[4];
        $time[5] += 1900;
    }
    else
    {
        $class->warning_message("failed to parse $date with strptime");
        # fall back to old method
        @time = reverse(split(m,[-\s:/],, $date));
    }
    return @time;
}

=pod

=item numbers_to_datetime

  $datetime = UR::Time->numbers_to_datetime($sec, $min, $hour
                                             $day, $mon, $year);

This method converts a list of numeric date and time data into a
date/time in the the default format.  In effect, this method is the
reverse of C<datetime_to_numbers> (see L<"datetime_to_numbers">).  If
an error occurs, C<undef> is returned.

=cut

# combine a bunch of date and time fields into a date in the default format
sub numbers_to_datetime
{
    my $class = shift;
    my ($s, $min, $h, $d, $mon, $y) = @_;

    # get date/time format
    my $format = $class->config('datetime');
    if ($format)
    {
        $class->debug_message("got date format %s", $format);
    }
    else
    {
        $class->error_message("no datetime format specified");
        return;
    }

    # return date time
    return POSIX::strftime($format, $s, $min, $h, $d, $mon - 1, $y - 1900);
}

=pod

=item datetime_to_time

  $time = UR::Time->datetime_to_time($datetime);

This method converts a date/time in the default format into epoch
seconds.

=cut

sub datetime_to_time
{
    my ($class, $date) = @_;

    # try to parse
    my @time = strptime($date);
    if (@time)
    {
        $class->debug_message("parse date $date");
    }
    else
    {
        $class->warning_message("failed to parse $date with strptime");
        # fall back to old method
        @time = reverse(split(m,[-\s:/],, $date));
        # fix month
        --$time[4];
        # fix year
        $time[5] -= 1900;
    }

    my $time = timelocal(@time);

    return $time;
}

=pod

=item add_date_delta_days
  
    $date = UR::Time->add_date_delta_days($date, $days);

This method accepts a date in the default format and will return a
date $days in the past or future.  Use negative values for calculating
dates in the past and positive values for dates in the future.  It is
based on C<Date::Pcalc::Add_Delta_Days> (see
L<Date::Pcalc/"Add_Delta_Days">).

=cut

sub add_date_delta_days 
{
    my $class = shift;
    my ($date, $days) = @_;

    # get date format
    my $format = $class->config('date');
    if ($format)
    {
        $class->debug_message("got date format %s", $format);
    }
    else
    {
        $class->error_message("no date format specified");
        return;
    }

    # try to parse date
    my ($d, $m, $y);
    my @time = strptime($date);
    if (@time)
    {
        $class->debug_message("parsed date $date: " . join('-', @time[3 .. 5]));
        ($d, $m, $y) = @time[3 .. 5];
        # correct for use with Date::Pcalc
        ++$m;
        $y += 1900;
    }
    else
    {
        $class->warning_message("failed to parse date with strptime");
        # fall back to old method
        ($y, $m, $d) = split(m/-/, $date);
    }

    # add the days
    ($y, $m, $d) = Date::Pcalc::Add_Delta_Days($y, $m, $d, $days);
    
    # format the string
    return POSIX::strftime($format, 0, 0, 0, $d, --$m, $y - 1900);
}

1;
__END__

=pod

=back

=head1 SEE ALSO

Date::Pcalc(3), Date::Parse(3), strftime(3), POSIX(3)

=cut

#$Header$
