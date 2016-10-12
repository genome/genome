package Genome::Command::Viewer;

use strict;
use warnings;

use Term::ReadKey 'GetTerminalSize';
use IO::Handle;

use Genome::Utility::Text qw(justify);

class Genome::Command::Viewer {
    doc => "Base class for 'view' commands.",
    is => 'Command::V2',
    has => [
        COLUMN_WIDTH => {
            is => 'Number',
            default_value => 47,
        }
    ],
};


sub write_report {
    my ($self, $width, $handle) = @_;
    Carp::croak("Method write_report must be implemented in sub-classes.");
}

sub get_terminal_width {
    # GetTerminalSize() returns an empty array if unsupported, e.g. in apipe-ci logs.
    local $SIG{__WARN__} = sub { };
    my ($screen_width) = GetTerminalSize();
    $screen_width ||= 80;
    return $screen_width;
}

sub execute {
    my ($self) = @_;


    my $handle = new IO::Handle;
    STDOUT->autoflush(1);
    $handle->fdopen(fileno(STDOUT), 'w');

    my $screen_width = $self->get_terminal_width();
    $self->write_report($screen_width, $handle);
    return 1;
}

sub get_report {
    my ($self, $width) = @_;

    my $handle = new IO::String;
    $self->write_report($width, $handle);
    my $report = ${$handle->string_ref};
    $handle->close();

    return $report;
}

sub _display_many {
    my ($self, $handle, $section_name, $method_name, @items) = @_;

    print $handle $self->_color_heading($section_name) . "\n";
    if(@items) {
        for my $item (@items) {
            $self->$method_name($handle, $item);
        }
    } else {
        print $handle "None\n";
    }
    print $handle "\n";
}

sub _software_result_test_names {
    my ($self, @results) = @_;

    unless (scalar(@results)) {
        return "No results found";
    }

    my $UNDEF = '***undef***';
    my @test_names = map {$_->test_name ? $_->test_name : $UNDEF} @results;

    my %counts;
    $counts{$_}++ for @test_names;
    my @sorted_counts = (
        sort { $b->[1] <=> $a->[1] }
        map { [ $_, $counts{$_} ] } keys %counts
    );

    my @entries;
    while (my $entry = shift(@sorted_counts)) {
        my ($name, $count) = @{$entry};
        if ($name eq $UNDEF) {
            $name = "undef";
        } else {
            $name = '"' . $name . '"';
        }

        push @entries, sprintf("%s (%d)", $name, $count);
    }
    return join(", ", @entries);
}

sub _clean_up_timestamp {
    my ($self, $dirty_stamp) = @_;

    unless (defined $dirty_stamp) {
        return '';
    }

    my ($clean_stamp) = split('\.', $dirty_stamp);
    $clean_stamp =~ s/\s$//;
    return $clean_stamp;
}

sub _pad_left {
    my ($self, $arg, $length) = @_;
    return justify($arg, 'right', $length);

}

sub _pad_right {
    my ($self, $arg, $length) = @_;
    return justify($arg, 'left', $length);
}

1;
