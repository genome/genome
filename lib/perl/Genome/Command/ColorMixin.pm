package Genome::Command::ColorMixin;

use strict;
use warnings;

class Genome::Command::ColorMixin {
    doc => "Mixin that gives commands color option",
    has => [
        color => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 1,
            doc => 'Display report in color.'
        },
    ],
};

sub _color {
    my $self = shift;
    my $string = shift;

    if(-t STDOUT and -t STDERR and $self->color and @_) {
        return Term::ANSIColor::colored($string, @_);
    } else {
        return $string;
    }
}

our %STATUS_COLORS = (
    new => "white",
    scheduled => "white",

    running => "cyan",
    'running*' => 'cyan',

    done      => "green",
    succeeded => "green",

    abandoned => "magenta",

    crashed     => "red",
    failed      => "red",
    unstartable => "red",
);

sub _status_colors {
    my $self = shift;
    my $status = shift;
    return $STATUS_COLORS{$status};
}

sub _status_color {
    my ($self, $text) = @_;
    return $self->_colorize_text_by_map($text, $text, %STATUS_COLORS);
}

sub _colorize_text_by_map {
    my ($self, $text, $color_key, %color_map) = @_;

    my $stripped_key = $self->_strip_key($color_key);
    if (exists $color_map{$stripped_key}) {
        return $self->_color($text, $color_map{$stripped_key});
    }

    return $text;
}

sub _strip_key {
    my ($self, $text) = @_;

    my $stripped_text = $text;
    $stripped_text =~ tr/A-Z/a-z/;
    $stripped_text =~ s/ //g;

    return $stripped_text;
}



1;
